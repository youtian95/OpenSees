// ============================================================================
// 文件: MasonryBendingMat.h
// 功能: MasonryBendingMat 砌体弯曲单轴材料类声明
// 说明: 由原 MasonryMat 拆分而来(剪切版本见 MasonryShearMat), 已完成注册:
//       - classTags.h: MAT_TAG_MasonryBendingMat (237)
//       - CMakeLists.txt / uniaxial/Makefile / SRC/Makefile: 构建系统
//       - OpenSeesUniaxialMaterialCommands.cpp: 命令分发(Tcl 与 OpenSeesPy 共用)
//       - TclModelBuilderUniaxialMaterialCommand.cpp: 经典 Tcl 分发
// 用法: 脚本命令
//       uniaxialMaterial MasonryBendingMat matTag? Ke? Mmax? <uu?> <R_My?> <alpha1?> <alpha2?> <gamma1?> <gamma2?>
//       Ke: 弹性刚度; Mmax: 最大弯矩; uu: 极限弯曲角度(论文取0.008rad, 默认0.008)
//       R_My: 屈服弯矩与Mmax的比值(论文取0.7)
//       alpha1: 第一段卸载段的力下降与下降点的力之比(论文取0.2)
//       alpha2: 第二段卸载段的力下降与下降点(alpha1*M_A)的力之比(论文取0.1)
//       gamma1: 第一段卸载段的刚度与Ke之比(论文取1.2)
//       gamma2: 第二段卸载段的刚度与屈服后刚度之比(论文取1.2)
//       实现见同目录 MasonryBendingMat.cpp
// ============================================================================

#ifndef MasonryBendingMat_h
#define MasonryBendingMat_h

#include <UniaxialMaterial.h>

// MasonryBendingMat: 砌体弯曲单轴材料
class MasonryBendingMat : public UniaxialMaterial
{
  public:
    // 构造函数
    // tag: 材料标签; Ke Mmax 必选; uu~gamma2 可选(默认值见参数声明)
    MasonryBendingMat(int tag, double Ke, double Mmax,
                      double uu = 0.008, double R_My = 0.7,
                      double alpha1 = 0.2, double alpha2 = 0.1,
                      double gamma1 = 1.2, double gamma2 = 1.2);
    MasonryBendingMat();   // 默认构造(FEM_ObjectBroker 并行/数据库恢复时需要)
    ~MasonryBendingMat();  // 析构函数

    // ---- UniaxialMaterial 纯虚函数(必须实现) ----
    // 设置试应变(与应变率), 在其中由本构规律计算应力与切线刚度
    int setTrialStrain(double strain, double strainRate = 0.0);
    double getStrain(void);          // 返回当前试应变
    double getStress(void);          // 返回当前试应力
    double getTangent(void);         // 返回当前试切线刚度
    double getInitialTangent(void);  // 返回初始刚度(Ke)

    int commitState(void);           // 提交状态: 试状态 -> 已提交状态
    int revertToLastCommit(void);    // 回滚: 用已提交状态恢复试状态
    int revertToStart(void);         // 回滚到初始状态
    UniaxialMaterial *getCopy(void); // 返回材料的完整副本

    // ---- 并行/数据库接口(MovableObject/TaggedObject 要求) ----
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);  // 打印材料信息

  private:
    // ---- 材料参数 ----
    double Ke;      // 弹性刚度
    double Mmax;    // 最大弯矩
    double uu;      // 极限弯曲角度
    double R_My;    // 屈服弯矩与Mmax的比值
    double alpha1;  // 第一段卸载段力下降比
    double alpha2;  // 第二段卸载段力下降比
    double gamma1;  // 第一段卸载段刚度与Ke之比
    double gamma2;  // 第二段卸载段刚度与屈服后刚度之比

    // ---- 试状态(当前迭代步) ----
    double trialStrain, trialStrainRate;  // 试应变与试应变率
    double trialStress, trialTangent;     // 试应力与试切线刚度

    // ---- 已提交状态(上一收敛步) ----
    double cStrain, cStress, cTangent;

    // TODO: 在此添加弯曲本构所需历史变量(如卸载/再加载标志、塑性转角等),
    //       并同步在 commitState/revertToStart/sendSelf/recvSelf 中处理
};

#endif