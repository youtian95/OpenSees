// ============================================================================
// 文件: MasonryBendingMat.h
// 功能: MasonryBendingMat 砌体弯曲单轴材料类声明
// 用法: 脚本命令
//       uniaxialMaterial MasonryBendingMat matTag? Ke? Mmax? <uu?> <R_My?> <CF?> <CD?> <gamma1?> <gamma2?> <postUltimateStiffness?>
//       Ke: 弹性刚度; Mmax: 最大弯矩; uu: 极限弯曲角度(论文取0.008rad, 默认0.008)
//       R_My: 屈服弯矩与Mmax的比值(论文取0.7)
//       CF: 第一段卸载相对反转点弯矩的下降比例(论文取0.2)
//       CD: 卸载第二折点相对屈服转角的水平位置参数(论文取0.1)
//       gamma1: 极限转角处第一段卸载刚度与Ke之比
//       gamma2: 第二段卸载路径斜率系数(论文取1.2，对应NextFEM的cC作用)
//       postUltimateStiffness: 超过极限转角后的骨架刚度，默认0；负值使承载力线性下降至0后保持为0
// 说明: 当前弯曲模型不包含剪切材料中的beta退化参数，不考虑强度或能量退化。
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
    // tag: 材料标签; Ke Mmax 必选; uu~postUltimateStiffness 可选(默认值见参数声明)
    MasonryBendingMat(int tag, double Ke, double Mmax, double uu = 0.008, double R_My = 0.7, double CF = 0.2, double CD = 0.1, double gamma1 = 1.2, double gamma2 = 1.2, double postUltimateStiffness = 0.0);
    MasonryBendingMat();   // 默认构造(FEM_ObjectBroker 并行/数据库恢复时需要)
    ~MasonryBendingMat();  // 析构函数

    // ---- UniaxialMaterial 纯虚函数(必须实现) ----
    // 设置试应变(与应变率), 在其中由本构规律计算应力与切线刚度
    int setTrialStrain(double strain, double strainRate = 0.0);
    // 以初始刚度计算线性试算状态，不推进弯曲材料的非线性历史。
    // strain: 当前端部转角。
    // strainRate: 当前端部转角速率。
    int setTrialLinearStrain(double strain, double strainRate = 0.0);
    double getStrain(void);          // 返回当前试应变
    double getStress(void);          // 返回当前试应力
    double getTangent(void);         // 返回当前试切线刚度
    double getInitialTangent(void);  // 返回初始刚度(Ke)
    double getYieldStrain(void) const; // 返回当前骨架的屈服转角。

    // 设置当前试算轴力对应的弯曲骨架峰值；Mmax: 当前试算最大弯矩
    int setTrialBackbone(double Mmax);

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
    double initialMmax;    // 创建材料时输入的初始最大弯矩
    double Mmax;           // 当前试算最大弯矩
    double committedMmax;  // 上一收敛状态的最大弯矩
    double uu;      // 极限弯曲角度
    double R_My;    // 屈服弯矩与Mmax的比值
    double CF;      // 第一段卸载段相对反转点的力下降比
    double CD;      // 第二折点相对屈服转角的水平位置参数
    double gamma1;  // 极限转角处第一段卸载刚度与Ke之比
    double gamma2;  // 第二段卸载路径斜率系数
    double postUltimateStiffness;  // 超过极限转角后的骨架刚度，负刚度下降至零承载力后截断

    // ---- 分支状态机 ----
    enum Branch : int {
      ELASTIC = 0,                         // 弹性段
      HARDENING_1 = 1,                     // 屈服后骨架段
      UNLOADING_1 = 5,                     // 第一段卸载
      UNLOADING_2 = 6,                     // 第二段卸载
      UNLOADING_3 = 7,                     // 第三段卸载至反向屈服点
      RELOADING_FROM_UNLOADING = 8         // 卸载后反向再加载
    };

    // ---- 试状态与已提交状态 ----
    // 所有状态变量均为值类型，保证默认拷贝能够复制完整材料历史。
    struct State {
      double strain = 0.0;                // 当前转角
      double strainRate = 0.0;            // 当前转角速度
      double stress = 0.0;                // 当前弯矩
      double tangent = 0.0;               // 当前切线刚度
      Branch branch = ELASTIC;            // 当前分支
      double ldir = 1.0;                  // 当前加载方向
      double revStrain = 0.0;             // 反转点A的转角
      double revStress = 0.0;             // 反转点A的弯矩
      double maxAbsStrain = 0.0;          // 已达到的历史最大绝对转角
      bool directUnloading3 = false;      // 是否从再加载直接进入第三段
    };

    State tState;                         // 试状态
    State cState;                         // 已提交状态

    // ---- 推算参数 ----
    double My;                             // 屈服弯矩
    double uy;                             // 屈服转角
    double Kp;                             // 屈服后的骨架刚度

    // ---- 辅助函数 ----
    // 根据当前试算最大弯矩重算屈服弯矩、屈服转角和屈服后刚度
    void updateDerivedParameters();
    // 根据历史最大绝对转角计算第一段卸载刚度，超过极限转角后继续线性退化。
    // state: 保存历史最大绝对转角的材料状态。
    double unloadingStiffness(const State &state) const;
    void backbone(double rotation, double &moment, double &tangent);
    State initialState() const;
    // 根据反转点弯矩方向计算第三段卸载的反向屈服转角。
    // state: 当前卸载状态，包含反转点弯矩和加载方向。
    double unloadingTargetStrain(const State &state) const;
    void computeUnloadingPoints(
        const State &state, double &unload1Strain, double &unload1Stress,
        double &unload2Strain, double &unload2Stress);
    void evaluateUnloadingPath(State &state, double rotation);

    bool checkTransitions();

    void ruleFromElastic();
    void ruleFromHardening1();
    void ruleFromUnloading1();
    void ruleFromUnloading2();
    void ruleFromUnloading3();
    void ruleFromReloading();
};

#endif
