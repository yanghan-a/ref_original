#include "6dof_kinematic.h"
#include "communication.hpp"
inline float cosf(float x)
{
    return arm_cos_f32(x);
}

inline float sinf(float x)
{
    return arm_sin_f32(x);
}

static void MatMultiply(const float* _matrix1, const float* _matrix2, float* _matrixOut,
                        const int _m, const int _l, const int _n)
{
    float tmp;
    int i, j, k;
    for (i = 0; i < _m; i++)
    {
        for (j = 0; j < _n; j++)
        {
            tmp = 0.0f;
            for (k = 0; k < _l; k++)
            {
                tmp += _matrix1[_l * i + k] * _matrix2[_n * k + j];
            }
            _matrixOut[_n * i + j] = tmp;
        }
    }
}

static void RotMatToEulerAngle(const float* _rotationM, float* _eulerAngles)
{//旋转矩阵转XYZ固定角
    float A, B, C, cb;

    if (fabs(_rotationM[6]) >= 1.0 - 0.0001)
    {
        if (_rotationM[6] < 0)
        {
            A = 0.0f;
            B = (float) M_PI_2;
            C = atan2f(_rotationM[1], _rotationM[4]);
        } else
        {
            A = 0.0f;
            B = -(float) M_PI_2;
            C = -atan2f(_rotationM[1], _rotationM[4]);
        }
    } else
    {
        B = atan2f(-_rotationM[6], sqrtf(_rotationM[0] * _rotationM[0] + _rotationM[3] * _rotationM[3]));//beta
        cb = cosf(B);
        A = atan2f(_rotationM[3] / cb, _rotationM[0] / cb);//alpha
        C = atan2f(_rotationM[7] / cb, _rotationM[8] / cb);//gama
    }

    _eulerAngles[0] = C;
    _eulerAngles[1] = B;
    _eulerAngles[2] = A;
}

static void EulerAngleToRotMat(const float* _eulerAngles, float* _rotationM)
{
    float ca, cb, cc, sa, sb, sc;

    cc = cosf(_eulerAngles[0]);
    cb = cosf(_eulerAngles[1]);
    ca = cosf(_eulerAngles[2]);
    sc = sinf(_eulerAngles[0]);
    sb = sinf(_eulerAngles[1]);
    sa = sinf(_eulerAngles[2]);

    _rotationM[0] = ca * cb;
    _rotationM[1] = ca * sb * sc - sa * cc;
    _rotationM[2] = ca * sb * cc + sa * sc;
    _rotationM[3] = sa * cb;
    _rotationM[4] = sa * sb * sc + ca * cc;
    _rotationM[5] = sa * sb * cc - ca * sc;
    _rotationM[6] = -sb;
    _rotationM[7] = cb * sc;
    _rotationM[8] = cb * cc;
}


DOF6Kinematic::DOF6Kinematic(float L_BS, float D_BS, float L_AM, float L_FA, float D_EW, float L_WT)
    : armConfig(ArmConfig_t{L_BS, D_BS, L_AM, L_FA, D_EW, L_WT})
{
    float tmp_DH_matrix[6][4] = {// theta, d, a, alpha
        {0.0f,            armConfig.L_BASE,    armConfig.D_BASE, -(float) M_PI_2},
        {-(float) M_PI_2, 0.0f,                armConfig.L_ARM,  0.0f},
        {0.0f,            0.0f,                armConfig.D_ELBOW, -(float) M_PI_2},//-->{0.0f,  0.0f,   armConfig.D_ELBOW, -(float) M_PI_2}
        {0.0f,            armConfig.L_FOREARM, 0.0f,             -(float) M_PI_2},
        {0.0f,            0.0f,                0.0f,             (float) M_PI_2},
        {0.0f,            armConfig.L_WRIST  , 0.0f,               0.0f}
    };
    memcpy(DH_matrix, tmp_DH_matrix, sizeof(tmp_DH_matrix));

    float tmp_L1_bs[3] = {armConfig.D_BASE, -armConfig.L_BASE, 0.0f};//这里的定义与一个关系有关就是坐标系A在B中的位置和B在A中的位置的关系
    memcpy(L1_base, tmp_L1_bs, sizeof(tmp_L1_bs));
    float tmp_L2_se[3] = {armConfig.L_ARM, 0.0f, 0.0f};
    memcpy(L2_arm, tmp_L2_se, sizeof(tmp_L2_se));
    float tmp_L3_ew[3] = {armConfig.D_ELBOW, 0.0f, 0.0f};
    memcpy(L3_elbow, tmp_L3_ew, sizeof(tmp_L3_ew));
    float tmp_L4_forearm[3] = {0.0f, -armConfig.L_FOREARM, 0.0f};
    memcpy(L4_forearm, tmp_L4_forearm, sizeof(tmp_L4_forearm));
    float tmp_L6_wt[3] = {0.0f, 0.0f, armConfig.L_WRIST};//注意这里加了负号
    memcpy(L6_wrist, tmp_L6_wt, sizeof(tmp_L6_wt));

    l_se_2 = armConfig.L_ARM * armConfig.L_ARM;
    l_se = armConfig.L_ARM;
    l_ew_2 = armConfig.L_FOREARM * armConfig.L_FOREARM + armConfig.D_ELBOW * armConfig.D_ELBOW;
    l_ew = sqrtf(l_ew_2);
    atan_e = atanf(armConfig.D_ELBOW / armConfig.L_FOREARM);
}

bool
DOF6Kinematic::SolveFK(const DOF6Kinematic::Joint6D_t &_inputJoint6D, DOF6Kinematic::Pose6D_t &_outputPose6D)
{
    float q_in[6];
    float q[6];
    float cosq, sinq;
    float cosa, sina;
    float P06[6];
    float R06[9];
    float R[6][9];//6个关节的旋转矩阵
    float R02[9];
    float R03[9];
    float R04[9];
    float R05[9];
    float L0_bs[3];
    float L0_se[3];
    float L0_ew[3];
    float L0_fa[3];
    float L0_wt[3];

    for (int i = 0; i < 6; i++)
        q_in[i] = _inputJoint6D.a[i] / RAD_TO_DEG;//把输入的关节角度转为弧度

    for (int i = 0; i < 6; i++)
    {
        q[i] = q_in[i];//q_in[i] + DH_matrix[i][0];
        cosq = cosf(q[i]);
        sinq = sinf(q[i]);
        cosa = cosf(DH_matrix[i][3]);
        sina = sinf(DH_matrix[i][3]);

        R[i][0] = cosq;
        R[i][1] = -cosa * sinq;
        R[i][2] = sina * sinq;
        R[i][3] = sinq;
        R[i][4] = cosa * cosq;
        R[i][5] = -sina * cosq;
        R[i][6] = 0.0f;
        R[i][7] = sina;
        R[i][8] = cosa;
    }

    MatMultiply(R[0], R[1], R02, 3, 3, 3);
    MatMultiply(R02, R[2], R03, 3, 3, 3);
    MatMultiply(R03, R[3], R04, 3, 3, 3);
    MatMultiply(R04, R[4], R05, 3, 3, 3);
    MatMultiply(R05, R[5], R06, 3, 3, 3);

    MatMultiply(R[0], L1_base, L0_bs, 3, 3, 1);
    MatMultiply(R02, L2_arm, L0_se, 3, 3, 1);
    MatMultiply(R03, L3_elbow, L0_ew, 3, 3, 1);
    MatMultiply(R04, L4_forearm, L0_fa, 3, 3, 1);
    MatMultiply(R06, L6_wrist, L0_wt, 3, 3, 1);

    for (int i = 0; i < 3; i++)
        P06[i] = L0_bs[i] + L0_se[i] + L0_ew[i] + L0_fa[i] + L0_wt[i];

    RotMatToEulerAngle(R06, &(P06[3]));

    _outputPose6D.X = P06[0];
    _outputPose6D.Y = P06[1];
    _outputPose6D.Z = P06[2];
    _outputPose6D.A = P06[3] * RAD_TO_DEG;
    _outputPose6D.B = P06[4] * RAD_TO_DEG;
    _outputPose6D.C = P06[5] * RAD_TO_DEG;
    memcpy(_outputPose6D.R, R06, 9 * sizeof(float));

    return true;
}

bool DOF6Kinematic::SolveIK(const DOF6Kinematic::Pose6D_t &_inputPose6D, const Joint6D_t &_lastJoint6D,
                            DOF6Kinematic::IKSolves_t &_outputSolves)
{
    float qs[2];//theta1
    float qa[2][2];//theta2、theta3
    float qw[2][3];//theta4、theta5、theta6
    float cosqs, sinqs;
    float cosqa[2], sinqa[2];
    float cosqw, sinqw;
    float P06[6];//六维姿态
    float R06[9];//旋转矩阵
    float P0_w[3];
    float P1_w[3];
    float L0_wt[3];
    float L1_sw[3];
    float R10[9];
    float R31[9];
    float R30[9];
    float R36[9];
    float l_sw_2, l_sw, atan_a, acos_a, acos_e;

    int ind_arm, ind_elbow, ind_wrist;
    int i;

    // if (0 == l_ew)
    // {
    //     l_ew = sqrtf(l_ew_2);
    //     atan_e = atanf(armConfig.D_ELBOW / armConfig.L_FOREARM);
    // }

    P06[0] = _inputPose6D.X / 1000.0f;//转换为m
    P06[1] = _inputPose6D.Y / 1000.0f;
    P06[2] = _inputPose6D.Z / 1000.0f;
    if (!_inputPose6D.hasR)//转换为弧度制
    {
        P06[3] = _inputPose6D.A / RAD_TO_DEG;
        P06[4] = _inputPose6D.B / RAD_TO_DEG;
        P06[5] = _inputPose6D.C / RAD_TO_DEG;
        EulerAngleToRotMat(&(P06[3]), R06);//计算XYZ固定角对应的旋转矩阵
    } else
    {
        memcpy(R06, _inputPose6D.R, 9 * sizeof(float));
    }
    for (i = 0; i < 2; i++)
    {
        qs[i] = _lastJoint6D.a[0];
        qa[i][0] = _lastJoint6D.a[1];
        qa[i][1] = _lastJoint6D.a[2];
        qw[i][0] = _lastJoint6D.a[3];
        qw[i][1] = _lastJoint6D.a[4];
        qw[i][2] = _lastJoint6D.a[5];
    }
    MatMultiply(R06, L6_wrist, L0_wt, 3, 3, 1);//计算最后一个L_WRIST在0坐标系中的向量
    for (i = 0; i < 3; i++)//计算关节4、5的坐标点
    {
        P0_w[i] = P06[i] - L0_wt[i];//用末端在0坐标系中的坐标减去link L_WRIST在0坐标系中的坐标，得到共轴点的位置
    }
    if (sqrt(P0_w[0] * P0_w[0] + P0_w[1] * P0_w[1]) <= 0.000001)//判断x、y坐标是否为0
    {//这个情况可以直接忽略，机械臂不可能到达这个位置
        qs[0] = _lastJoint6D.a[0];//此时theta1可以为任意值，为了减少运动，那么就认定位置不变
        qs[1] = _lastJoint6D.a[0];
        for (i = 0; i < 4; i++)
        {
            _outputSolves.solFlag[0 + i][0] = -1;//-1代表theta1可以为任意值，为了减少运动，那么就认定位置不变
            _outputSolves.solFlag[4 + i][0] = -1;
        }
    } else
    {
        qs[0] = atan2f(P0_w[1], P0_w[0]);//这里直接计算出theta1的角度
        qs[1] = atan2f(-P0_w[1], -P0_w[0]);//theta1的两个解，两个角相差180度，事实上这个解得舍去
        for (i = 0; i < 4; i++)
        {
            _outputSolves.solFlag[0 + i][0] = 1;//1代表正常解
            _outputSolves.solFlag[4 + i][0] = 1;
        }
    }
    for (ind_arm = 0; ind_arm < 2; ind_arm++)
    {//下面一部分是求解theta2、theta3，主要运用的就是几何法，对两个三角形考虑余弦定理
        cosqs = cosf(qs[ind_arm] );
        sinqs = sinf(qs[ind_arm]);

        R10[0] = cosqs;//构造R01矩阵的逆，也就是R10矩阵的转置
        R10[1] = sinqs;
        R10[2] = 0.0f;
        R10[3] = 0.0f;
        R10[4] = 0.0f;
        R10[5] = -1.0f;
        R10[6] = -sinqs;
        R10[7] = cosqs;
        R10[8] = 0.0f;

        MatMultiply(R10, P0_w, P1_w, 3, 3, 1);
        for (i = 0; i < 3; i++)//在坐标系1中计算
        {
            L1_sw[i] = P1_w[i] - L1_base[i];//这里要改为+，因为定义的时候添加了负号
        }
        l_sw_2 = L1_sw[0] * L1_sw[0] + L1_sw[1] * L1_sw[1];
        l_sw = sqrtf(l_sw_2);

        if (fabs(l_se + l_ew - l_sw) <= 0.000001)//这里是判断是否是边界
        {
            qa[0][0] = atan2f(L1_sw[1], L1_sw[0]);//theta2
            qa[1][0] = qa[0][0];
            qa[0][1] = -(((float) M_PI_2) - atan_e);
            qa[1][1] = -(((float) M_PI_2) - atan_e);
            if (l_sw > l_se + l_ew)
            {//超出范围
                for (i = 0; i < 2; i++)
                {
                    _outputSolves.solFlag[4 * ind_arm + 0 + i][1] = 0;//0代表到不了，但几乎重合
                    _outputSolves.solFlag[4 * ind_arm + 2 + i][1] = 0;
                }
            } else
            {
                for (i = 0; i < 2; i++)
                {
                    _outputSolves.solFlag[4 * ind_arm + 0 + i][1] = 1;//1代表正常解
                    _outputSolves.solFlag[4 * ind_arm + 2 + i][1] = 1;
                }
            }
        } else if (fabs(l_sw - fabs(l_se - l_ew)) <= 0.000001)//另一个边界，这个情况机械臂也达不到位置
        {
            qa[0][0] = atan2f(L1_sw[1], L1_sw[0]);
            qa[1][0] = qa[0][0];
            if (0 == ind_arm)
            {
                qa[0][1] = (float) M_PI_2+ atan_e;
                qa[1][1] = -(float) M_PI_2*3+ atan_e;
            } else
            {
                qa[0][1] = (float) M_PI_2+ atan_e;
                qa[1][1] = -(float) M_PI_2*3+ atan_e;
            }
            if (l_sw < fabs(l_se - l_ew))
            {
                for (i = 0; i < 2; i++)
                {
                    _outputSolves.solFlag[4 * ind_arm + 0 + i][1] = 0;
                    _outputSolves.solFlag[4 * ind_arm + 2 + i][1] = 0;
                }
            } else
            {
                for (i = 0; i < 2; i++)
                {
                    _outputSolves.solFlag[4 * ind_arm + 0 + i][1] = 1;
                    _outputSolves.solFlag[4 * ind_arm + 2 + i][1] = 1;
                }
            }
        } else
        {
            atan_a = atan2f(L1_sw[1], L1_sw[0]);
            acos_a = 0.5f * (l_se_2 + l_sw_2 - l_ew_2) / (l_se * l_sw);
            if (acos_a >= 1.0f) acos_a = 0.0f;//即使超出范围，它也会保证关节二四的连线与4的目标点共线
            else if (acos_a <= -1.0f) acos_a = 0.0f;//二次方程角度考虑，这个情况不可能
            else acos_a = acosf(acos_a);
            acos_e = 0.5f * (l_se_2 + l_ew_2 - l_sw_2) / (l_se * l_ew);
            if (acos_e >= 1.0f) acos_e = 0.0f;
            else if (acos_e <= -1.0f) acos_e = (float) M_PI;
            else acos_e = acosf(acos_e);
            if (0 == ind_arm)
            {
                qa[0][0] = atan_a - acos_a ;
                qa[0][1] = atan_e - acos_e + (float) M_PI_2;
                qa[1][0] = atan_a + acos_a ;
                qa[1][1] = atan_e + acos_e - (float) M_PI_2*3;

            } else
            {
                qa[0][0] = atan_a + acos_a ;
                qa[0][1] = -atan_e + acos_e - (float) M_PI_2;
                qa[1][0] = atan_a - acos_a ;
                qa[1][1] = -atan_e - acos_e + (float) M_PI_2*3;
            }
            for (i = 0; i < 2; i++)
            {
                _outputSolves.solFlag[4 * ind_arm + 0 + i][1] = 1;
                _outputSolves.solFlag[4 * ind_arm + 2 + i][1] = 1;
            }
        }

        //下面开始算后三个关节，与姿态相关
        for (ind_elbow = 0; ind_elbow < 2; ind_elbow++)
        {

            cosqa[0] = cosf(qa[ind_elbow][0] );
            sinqa[0] = sinf(qa[ind_elbow][0] );
            cosqa[1] = cosf(qa[ind_elbow][1] );
            sinqa[1] = sinf(qa[ind_elbow][1] );

            R31[0] = cosqa[0] * cosqa[1] - sinqa[0] * sinqa[1];//这里根据已经求出来的theta2和theta3，计算R31矩阵，也就是R13的逆
            R31[1] = cosqa[0] * sinqa[1] + sinqa[0] * cosqa[1];
            R31[2] = 0.0f;
            R31[3] = 0.0f;
            R31[4] = 0.0f;
            R31[5] = -1.0f;
            R31[6] = -cosqa[0] * sinqa[1] - sinqa[0] * cosqa[1];
            R31[7] = cosqa[0] * cosqa[1] - sinqa[0] * sinqa[1];
            R31[8] = 0.0f;

            MatMultiply(R31, R10, R30, 3, 3, 3);
            MatMultiply(R30, R06, R36, 3, 3, 3);//计算R36矩阵，然后ZYZ角解出theta4、5、6

            // printf("R36:,%f,%f,%f,%f,%f,%f,%f,%f,%f\r\n", R36[0],R36[1],R36[2],R36[3],R36[4],R36[5],R36[6],R36[7],R36[8]);
            if (R36[8] >= 1.0 - 0.000001)
            {
                cosqw = 1.0f;
                qw[0][1] = 0.0f;
                qw[1][1] = 0.0f;
            } else if (R36[8] <= -1.0 + 0.000001)
            {
                cosqw = -1.0f;
                if (0 == ind_arm)
                {
                    qw[0][1] = (float) M_PI;
                    qw[1][1] = -(float) M_PI;
                } else
                {
                    qw[0][1] = -(float) M_PI;
                    qw[1][1] = (float) M_PI;
                }
            } else
            {
                cosqw = R36[8];
                if (0 == ind_arm)
                {
                    qw[0][1] = acosf(cosqw);
                    qw[1][1] = -acosf(cosqw);
                } else
                {
                    qw[0][1] = -acosf(cosqw);
                    qw[1][1] = acosf(cosqw);
                }
            }
            if (1.0f == cosqw || -1.0f == cosqw)
            {
                if (0 == ind_arm)
                {
                    qw[0][0] = atan2f(R36[3], R36[0]);
                    qw[0][2] = 0.0f;
                    qw[1][2] = 0.0f;
                    qw[1][0] = qw[0][0];
                } else
                {
                    qw[0][0] = atan2f(R36[3], R36[0]);
                    qw[0][2] = 0.0f;
                    qw[1][2] = 0.0f;
                    qw[1][0] = qw[0][0];
                }
                _outputSolves.solFlag[4 * ind_arm + 2 * ind_elbow + 0][2] = -1;
                _outputSolves.solFlag[4 * ind_arm + 2 * ind_elbow + 1][2] = -1;
            } else
            {
                if (0 == ind_arm)
                {
                    qw[0][0] = atan2f(R36[5], R36[2]);
                    qw[1][0] = atan2f(-R36[5], -R36[2]);
                    qw[0][2] = atan2f(R36[7], -R36[6]);
                    qw[1][2] = atan2f(-R36[7], R36[6]);
                } else
                {
                    qw[0][0] = atan2f(-R36[5], -R36[2]);
                    qw[1][0] = atan2f(R36[5], R36[2]);
                    qw[0][2] = atan2f(-R36[7], R36[6]);
                    qw[1][2] = atan2f(R36[7], -R36[6]);
                }
                _outputSolves.solFlag[4 * ind_arm + 2 * ind_elbow + 0][2] = 1;
                _outputSolves.solFlag[4 * ind_arm + 2 * ind_elbow + 1][2] = 1;
            }
            for (ind_wrist = 0; ind_wrist < 2; ind_wrist++)
            {
                if (qs[ind_arm] > (float) M_PI)
                    _outputSolves.config[4 * ind_arm + 2 * ind_elbow + ind_wrist].a[0] =
                        qs[ind_arm] - (float) M_PI;
                else if (qs[ind_arm] < -(float) M_PI)
                    _outputSolves.config[4 * ind_arm + 2 * ind_elbow + ind_wrist].a[0] =
                        qs[ind_arm] + (float) M_PI;
                else
                    _outputSolves.config[4 * ind_arm + 2 * ind_elbow + ind_wrist].a[0] = qs[ind_arm];

                for (i = 0; i < 2; i++)
                {
                    if (qa[ind_elbow][i] > (float) M_PI)
                        _outputSolves.config[4 * ind_arm + 2 * ind_elbow + ind_wrist].a[1 + i] =
                            qa[ind_elbow][i] - (float) M_PI;
                    else if (qa[ind_elbow][i] < -(float) M_PI)
                        _outputSolves.config[4 * ind_arm + 2 * ind_elbow + ind_wrist].a[1 + i] =
                            qa[ind_elbow][i] + (float) M_PI;
                    else
                        _outputSolves.config[4 * ind_arm + 2 * ind_elbow + ind_wrist].a[1 + i] =
                            qa[ind_elbow][i];
                }

                for (i = 0; i < 3; i++)
                {
                    if (qw[ind_wrist][i] > (float) M_PI)
                        _outputSolves.config[4 * ind_arm + 2 * ind_elbow + ind_wrist].a[3 + i] =
                            qw[ind_wrist][i] - (float) M_PI;
                    else if (qw[ind_wrist][i] < -(float) M_PI)
                        _outputSolves.config[4 * ind_arm + 2 * ind_elbow + ind_wrist].a[3 + i] =
                            qw[ind_wrist][i] + (float) M_PI;
                    else
                        _outputSolves.config[4 * ind_arm + 2 * ind_elbow + ind_wrist].a[3 + i] =
                            qw[ind_wrist][i];
                }
            }
        }
    }

    for (i = 0; i < 8; i++)
        for (float &j: _outputSolves.config[i].a)
            j *= RAD_TO_DEG;

    return true;
}

DOF6Kinematic::Joint6D_t
operator-(const DOF6Kinematic::Joint6D_t &_joints1, const DOF6Kinematic::Joint6D_t &_joints2)
{
    DOF6Kinematic::Joint6D_t tmp{};
    for (int i = 0; i < 6; i++)
        tmp.a[i] = _joints1.a[i] - _joints2.a[i];

    return tmp;
}
