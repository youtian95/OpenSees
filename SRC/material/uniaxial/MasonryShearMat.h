// ============================================================================
// 文件: MasonryShearMat.h
// 功能: MasonryShearMat 砌体剪切单轴材料类声明
// 用法: 脚本命令
//       uniaxialMaterial MasonryShearMat matTag? Ke? Vmax? uu? <R_Vy?> <R_Vu?> <R_umax?> <alpha?> <gamma?> <beta?>
//       Ke: 弹性刚度; Vmax: 最大剪切力; uu: 极限剪切位移(论文取0.004h)
//       R_Vy: 屈服力与Vmax的比值(默认0.7)
//       R_Vu: 极限位移对应剪力与Vmax比值(默认0.8)
//       R_umax: 最大剪切力的位移与uu的比值(论文默认取0.5)
//       alpha: 卸载刚度除以弹性刚度(论文取0.8)
//       gamma: 卸载转折点的力与卸载点的力之比(默认取0.6, 论文中未出现该数值)
//       beta: 强度退化系数(论文取0.06), d_u2 = beta * d_E_h / Vmax
//       实现见同目录 MasonryShearMat.cpp
// ============================================================================

#ifndef MasonryShearMat_h
#define MasonryShearMat_h

#include <UniaxialMaterial.h>

// MasonryShearMat: 砌体剪切单轴材料
class MasonryShearMat : public UniaxialMaterial
{
  public:
    // 构造函数
    // tag: 材料标签; Ke Vmax uu 必选; R_Vy~beta 可选(默认值见参数声明)
    MasonryShearMat(int tag, double Ke, double Vmax, double uu, double R_Vy = 0.7, double R_Vu = 0.8, double R_umax = 0.5, double alpha = 0.8, double gamma = 0.6, double beta = 0.06);
    MasonryShearMat();   // 默认构造(FEM_ObjectBroker 并行/数据库恢复时需要)
    ~MasonryShearMat();  // 析构函数

    // ---- UniaxialMaterial 纯虚函数(必须实现) ----
    // 设置试应变(与应变率), 在其中由本构规律计算应力与切线刚度
    int setTrialStrain(double strain, double strainRate = 0.0);
    // 以初始刚度计算线性试算状态，不推进剪切材料的非线性历史。
    // strain: 当前剪切变形。
    // strainRate: 当前剪切变形速率。
    int setTrialLinearStrain(double strain, double strainRate = 0.0);
    double getStrain(void);          // 返回当前试应变
    double getStress(void);          // 返回当前试应力
    double getTangent(void);         // 返回当前试切线刚度
    double getInitialTangent(void);  // 返回初始刚度(Ke)
    double getYieldStrain(void) const; // 返回当前骨架的屈服位移。

    // 设置当前试算轴力对应的剪切骨架峰值；Vmax: 当前试算最大剪切力
    int setTrialBackbone(double Vmax);

    int commitState(void);           // 提交状态: 试状态 -> 已提交状态
    int revertToLastCommit(void);    // 回滚: 用已提交状态恢复试状态
    int revertToStart(void);         // 回滚到初始状态
    UniaxialMaterial *getCopy(void); // 返回材料的完整副本

    // ---- 并行/数据库接口(MovableObject/TaggedObject 要求) ----
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);  // 打印材料信息
    Response *setResponse(const char **argv, int argc, OPS_Stream &output);  // 输出材料内部诊断状态
    int getResponse(int responseID, Information &information);  // 返回材料内部诊断状态

  private:
    // ---- 材料参数 ----
    double Ke;      // 弹性刚度
    double initialVmax;    // 创建材料时输入的初始最大剪切力
    double Vmax;           // 当前试算最大剪切力
    double committedVmax;  // 上一收敛状态的最大剪切力
    double uu;      // 极限剪切位移
    double R_Vy;    // 屈服力与Vmax的比值
    double R_Vu;    // 极限位移对应剪力与Vmax比值
    double R_umax;  // 最大剪切力对应位移与uu的比值
    double alpha;   // 极限位移点卸载刚度与弹性刚度之比
    double gamma;   // 卸载转折点力与卸载点力之比
    double beta;    // 强度退化系数

    // ---- 推算的材料参数(可在构造函数中计算) ----
    double Vy;      // 屈服力
    double Vu;      // 极限位移对应剪力
    double umax;    // 最大剪切力对应位移
    double uy;     // 屈服位移

    // ---- 分支状态机: 穷举曲线中的每一条分支 ----
    enum Branch : int {
        ELASTIC = 0,   // 弹性阶段
        HARDENING_1 = 1,   // 硬化1阶段
        HARDENING_2 = 2,   // 硬化2阶段
        UNLOADING_1 = 5,     // 卸载1阶段
        UNLOADING_2 = 6,     // 卸载2阶段(反向加载)
        RELOADING_FROM_UNLOADING_1 = 7,  // 卸载1阶段再加载
    };

    // ---- 材料状态: 应变/应力/切线 + 全部历史变量 ----
    // 要求: 只放值类型(double/int), 不放指针, 保证默认拷贝安全
    struct State {
        double strain = 0.0;      // 应变
        double strainRate = 0.0;  // 应变率
        double stress = 0.0;      // 应力
        double tangent = 0.0;     // 切线刚度
        // 当前状态
        Branch branch = ELASTIC;  // 当前分支
        // 加载方向(这一步减去上一步的应变, 判断加载方向)
        double ldir = 1.0;     // 加载方向: 1.0 表示正向加载, -1.0 表示负向加载
        // 历史变量
        double umaxPos = 0.0;   // 历史最大绝对位移的正向值
        double umaxNeg = 0.0;   // 历史最大绝对位移的负向值，与umaxPos对称
        double revStrain = 0.0, revStress = 0.0;      // 反转点: 卸载线的起点
        double tgtStrain = 0.0, tgtStress = 0.0;  // 反向加载的目标点
        double cycleStartStrain = 0.0;  // 当前耗能回路起点位移
        double cycleStartStress = 0.0;  // 当前耗能回路起点剪力
        double cycleStartDir = 0.0;     // 当前耗能回路是否激活及其起点方向
        double cycleWork = 0.0;         // 从回路起点至当前状态的有符号外力功
    };

    State tState;    // 试状态(当前迭代步)
    State cState;   // 已提交状态(上一收敛步)

    // ---- 其他辅助私有函数 ----
    // 根据当前试算最大剪切力重算屈服力、残余力和特征位移
    void updateDerivedParameters();
    // 骨架曲线计算函数
    void backbone(double u, double &V, double &Et);
    // 计算射线与骨架曲线的交点
    void intersectionWithBackbone(double Strain, double Stress, double Tangent, double StrainDir, double &tgtStrain, double &tgtStress);
    // 计算卸载段1的刚度
    double computeUnload1Tangent(const State &state);
    // 计算从零点到指定骨架位移的功
    double backboneWork(double strain) const;
    // 根据原历史最大位移和当前回路耗能设置第二卸载阶段的退化目标
    void setUnloading2Target(State &state, const State &origin);
    // 将当前试增量的外力功累计到回路耗能状态
    void accumulateCycleWork();

    // 状态清零函数: 返回一个全零的状态, 仅切线刚度为 Ke
    State initialState() const;

    // 更新试状态是否变化
    bool checkTransitions();

    // 各个计算函数
    void ruleFromElastic(); 
    void ruleFromHardening1();
    void ruleFromHardening2();
    void ruleFromUnloading1();
    void ruleFromUnloading2();
    void ruleFromReloadingFromUnloading1();
};

#endif
