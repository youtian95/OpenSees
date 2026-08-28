// ============================================================================
// 文件: MasonryShearMat.cpp
// 功能: MasonryShearMat 砌体剪切单轴材料实现 + OPS 工厂函数
// 说明: 样板代码(工厂函数/状态管理/并行通讯)已完成,
//       核心本构规律在 setTrialStrain 中以 TODO 标出, 由开发者手动实现
// 用法: 编译后在脚本中使用:
//       Tcl:        uniaxialMaterial MasonryShearMat matTag? Ke? Vmax? uu? <可选参数> <-SymmetricMaxHistory>
//       OpenSeesPy: ops.uniaxialMaterial('MasonryShearMat', matTag, Ke, Vmax, uu, ...)
// 输入: matTag(整数); Ke Vmax uu(必选实数); R_Vy~beta(可选, 尾部省略用默认值); -SymmetricMaxHistory(可选)
// 输出: 创建 MasonryShearMat 材料对象并注册到模型构建器
// ============================================================================

#include "MasonryShearMat.h"

#include <OPS_Globals.h>
#include <elementAPI.h>
#include <classTags.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <Vector.h>
#include <algorithm>
#include <cmath>
#include <cstring>

// 已创建的材料实例数(用于只在首次创建时打印署名横幅)
static int numMasonryShearMat = 0;

// 工厂函数: uniaxialMaterial 命令的入口
// 命令语法(< > 内为可选参数, 尾部省略时使用构造函数的默认参数):
//   uniaxialMaterial MasonryShearMat matTag? Ke? Vmax? uu? <R_Vy?> <R_Vu?> <R_umax?> <alpha?> <gamma?> <beta?> <-SymmetricMaxHistory>
void *OPS_MasonryShearMat(void)
{
  // 首次创建时打印一次署名信息
  if (numMasonryShearMat == 0) {
    opserr << "MasonryShearMat uniaxial material - Written by youtian95\n";
    numMasonryShearMat = 1;
  }

  // 至少需要 matTag + Ke + Vmax + uu 四个参数
  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "WARNING invalid args, want: uniaxialMaterial MasonryShearMat tag? Ke? Vmax? uu? <R_Vy?> <R_Vu?> <R_umax?> <alpha?> <gamma?> <beta?> <-SymmetricMaxHistory>\n";
    return 0;
  }

  // 1) 读取整数参数: matTag
  int iData[1];
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial MasonryShearMat tag" << endln;
    return 0;
  }

  // 2) 依次读取位置数值参数和对称历史开关，未给出的尾部参数使用默认值。
  // OPS_GetStringFromAll 同时支持 Tcl 字符串和 OpenSeesPy 数值对象。
  double dData[9] = {0.0, 0.0, 0.0, 0.7, 0.8, 0.5, 0.8, 0.6, 0.06};
  int numberOfDoubleData = 0;
  bool symmetricMaxHistory = false;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    char argumentBuffer[128];
    const char *argument = OPS_GetStringFromAll(argumentBuffer, 128);
    if (std::strcmp(argument, "-SymmetricMaxHistory") == 0) {
      symmetricMaxHistory = true;
      continue;
    }
    if (numberOfDoubleData >= 9) {
      opserr << "WARNING too many MasonryShearMat numeric parameters\n";
      return 0;
    }
    char *end = 0;
    dData[numberOfDoubleData] = std::strtod(argument, &end);
    if (end == argument || *end != '\0') {
      opserr << "WARNING invalid MasonryShearMat parameter: " << argument << endln;
      return 0;
    }
    ++numberOfDoubleData;
  }
  if (numberOfDoubleData < 3) {
    opserr << "WARNING MasonryShearMat requires Ke, Vmax and uu\n";
    return 0;
  }

  // 3) 用完整参数和历史选项创建材料。
  UniaxialMaterial *theMaterial = new MasonryShearMat(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], dData[8], symmetricMaxHistory);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type MasonryShearMat\n";
    return 0;
  }
  return theMaterial;
}

// 全参构造函数
// tag: 材料标签; 其余参数含义见头文件
// 基类必须用 (tag, MAT_TAG_MasonryShearMat) 初始化
MasonryShearMat::MasonryShearMat(int tag, double _Ke, double _Vmax, double _uu, double _R_Vy, double _R_Vu, double _R_umax, double _alpha, double _gamma, double _beta, bool _symmetricMaxHistory)
  : UniaxialMaterial(tag, MAT_TAG_MasonryShearMat),
    Ke(_Ke), initialVmax(_Vmax), Vmax(_Vmax), committedVmax(_Vmax), uu(_uu),
    R_Vy(_R_Vy), R_Vu(_R_Vu), R_umax(_R_umax),
    alpha(_alpha), gamma(_gamma), beta(_beta), symmetricMaxHistory(_symmetricMaxHistory)
{
  // 根据初始骨架峰值计算屈服力、残余力和特征位移。
  updateDerivedParameters();
  // 历史变量初始化为零
  tState = initialState();
  cState = initialState();
}

// 默认构造函数: 供 FEM_ObjectBroker 在并行/数据库恢复时使用
MasonryShearMat::MasonryShearMat()
  : UniaxialMaterial(0, MAT_TAG_MasonryShearMat),
    Ke(0.0), initialVmax(0.0), Vmax(0.0), committedVmax(0.0), uu(0.0),
    R_Vy(0.7), R_Vu(0.8), R_umax(0.5),
    alpha(0.8), gamma(0.6), beta(0.06), symmetricMaxHistory(false),
    Vy(0.0), Vu(0.0), umax(0.0), uy(0.0)
{
  // 派生参数清零(recvSelf 恢复参数后会重算)
  Vy = 0.0; Vu = 0.0; umax = 0.0; uy = 0.0;
  // 状态初始化为零(切线刚度为 Ke, 此处 Ke=0)
  tState = initialState();
  cState = initialState();
}

// 析构函数
MasonryShearMat::~MasonryShearMat()
{
}

// 设置试应变并计算应力与切线刚度
// strain: 试应变; strainRate: 试应变率
// ★★★ 核心本构规律写在这里: 由 trialStrain 与历史变量计算 trialStress/trialTangent ★★★
int
MasonryShearMat::setTrialStrain(double strain, double strainRate)
{
  tState = cState;
  tState.strain = strain;
  tState.strainRate = strainRate;
  if (symmetricMaxHistory) {
    const double maximumAbsoluteStrain = std::max(std::max(cState.umaxPos, -cState.umaxNeg), std::abs(strain));
    tState.umaxPos = maximumAbsoluteStrain;
    tState.umaxNeg = -maximumAbsoluteStrain;
  } else {
    tState.umaxPos = std::max(cState.umaxPos, strain);
    tState.umaxNeg = std::min(cState.umaxNeg, strain);
  }

  // 计算加载方向
  double du = tState.strain - cState.strain;
  if (du > 1.0e-14)
      tState.ldir = 1.0;
  else if (du < -1.0e-14)
      tState.ldir = -1.0;
  else
      tState.ldir = cState.ldir;

  // 先根据当前试应变和历史状态确定下一步分支
  checkTransitions();

  // 各规则函数假设 tState.branch 已经是下一步状态
  switch (cState.branch) {
    case ELASTIC: ruleFromElastic(); break;
    case HARDENING_1: ruleFromHardening1(); break;
    case HARDENING_2: ruleFromHardening2(); break;
    case UNLOADING_1: ruleFromUnloading1(); break;
    case UNLOADING_2: ruleFromUnloading2(); break;
    case RELOADING_FROM_UNLOADING_1: ruleFromReloadingFromUnloading1(); break;
  }

  // 应力计算完成后再积分当前试增量的外力功，提交失败时会随试状态一起回滚。
  accumulateCycleWork();

  return 0;
}

int MasonryShearMat::setTrialLinearStrain(double strain, double strainRate)
{
  // 非控制弹簧始终沿原点处的初始弹性直线响应，同时保留上一步已提交历史。
  tState = cState;
  tState.strain = strain;
  tState.strainRate = strainRate;
  tState.stress = Ke * strain;
  tState.tangent = Ke;
  tState.branch = ELASTIC;
  return 0;
}

// 以下四个 get 函数返回当前试状态
double MasonryShearMat::getStrain(void)         { return tState.strain; }
double MasonryShearMat::getStress(void)         { return tState.stress; }
double MasonryShearMat::getTangent(void)        { return tState.tangent; }
double MasonryShearMat::getInitialTangent(void) { return Ke; }

double MasonryShearMat::getYieldStrain(void) const { return uy; }

// 设置当前试算轴力对应的剪切骨架峰值。
// Vmax: 当前试算最大剪切力；应在本次 setTrialStrain() 之前调用。
int MasonryShearMat::setTrialBackbone(double trialVmax)
{
  // 用初始承载力的百万分之一作为数值下限，避免零承载力导致派生参数除零。
  const double minimumVmax = 1.0e-6 * initialVmax;
  Vmax = std::max(trialVmax, minimumVmax);
  updateDerivedParameters();
  return 0;
}

// 提交状态: 试状态转为已提交状态
int
MasonryShearMat::commitState(void)
{
  cState = tState;
  committedVmax = Vmax;
  return 0;
}

// 回滚到上次提交状态(迭代不收敛时调用)
int
MasonryShearMat::revertToLastCommit(void)
{
  tState = cState;
  Vmax = committedVmax;
  updateDerivedParameters();
  return 0;
}

// 回滚到初始状态(重新开始分析时调用)
int
MasonryShearMat::revertToStart(void)
{
  Vmax = initialVmax;
  committedVmax = initialVmax;
  updateDerivedParameters();
  tState = initialState();
  cState = initialState();
  return 0;
}

// 返回材料完整副本(单元迭代时调用)
// 成员均为值类型, 默认拷贝构造即可; 若后续加入指针成员需自行深拷贝
UniaxialMaterial *
MasonryShearMat::getCopy(void)
{
  return new MasonryShearMat(*this);
}

// 打包发送(并行/数据库): 材料tag + 材料参数 + 已提交骨架峰值 + 完整已提交状态
// commitTag: 提交标签; theChannel: 通讯通道
int
MasonryShearMat::sendSelf(int commitTag, Channel &theChannel)
{
  // 打包布局: data(0)=材料tag; data(1..9)=材料参数; data(10)=已提交骨架峰值; data(11..26)=已提交状态; data(27)=对称历史开关
  static Vector data(28);
  data(0) = this->getTag();
  // 材料参数
  data(1) = Ke;     data(2) = initialVmax;   data(3) = uu;
  data(4) = R_Vy;   data(5) = R_Vu;   data(6) = R_umax;
  data(7) = alpha;  data(8) = gamma;  data(9) = beta;
  data(10) = committedVmax;
  // 已提交状态: 基本量 + 当前分支 + 全部历史变量
  data(11) = cState.strain;
  data(12) = cState.strainRate;
  data(13) = cState.stress;
  data(14) = cState.tangent;
  data(15) = (int)cState.branch;
  data(16) = cState.ldir;
  data(17) = cState.umaxPos;
  data(18) = cState.umaxNeg;
  data(19) = cState.revStrain;
  data(20) = cState.revStress;
  data(21) = cState.tgtStrain;
  data(22) = cState.tgtStress;
  data(23) = cState.cycleStartStrain;
  data(24) = cState.cycleStartStress;
  data(25) = cState.cycleStartDir;
  data(26) = cState.cycleWork;
  data(27) = symmetricMaxHistory ? 1.0 : 0.0;

  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0)
    opserr << "MasonryShearMat::sendSelf() - failed to send data\n";
  return res;
}

// 接收恢复(并行/数据库): 打包顺序必须与 sendSelf 严格对应
int
MasonryShearMat::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static Vector data(28);
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "MasonryShearMat::recvSelf() - failed to receive data\n";
    return res;
  }
  this->setTag((int)data(0));
  // 材料参数
  Ke = data(1);     initialVmax = data(2);   uu = data(3);
  R_Vy = data(4);   R_Vu = data(5);   R_umax = data(6);
  alpha = data(7);  gamma = data(8);  beta = data(9);
  committedVmax = data(10);
  Vmax = committedVmax;
  // 重算派生参数(派生量不传输, 必须由当前已提交骨架重建)
  updateDerivedParameters();
  // 已提交状态: 基本量 + 当前分支 + 全部历史变量
  cState.strain     = data(11);
  cState.strainRate = data(12);
  cState.stress     = data(13);
  cState.tangent    = data(14);
  cState.branch     = (Branch)(int)data(15);
  cState.ldir       = data(16);
  cState.umaxPos    = data(17);
  cState.umaxNeg    = data(18);
  cState.revStrain  = data(19);
  cState.revStress  = data(20);
  cState.tgtStrain  = data(21);
  cState.tgtStress  = data(22);
  cState.cycleStartStrain = data(23);
  cState.cycleStartStress = data(24);
  cState.cycleStartDir = data(25);
  cState.cycleWork = data(26);
  symmetricMaxHistory = data(27) != 0.0;
  // 用已提交状态初始化试状态
  tState = cState;
  return 0;
}

// 打印材料信息
// s: 输出流; flag: 输出详细程度标志
void
MasonryShearMat::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    // 按 OpenSees print -JSON 约定输出单个材料对象，布尔选项使用 JSON 布尔值。
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"MasonryShearMat\", ";
    s << "\"Ke\": " << Ke << ", ";
    s << "\"Vmax\": " << initialVmax << ", ";
    s << "\"uu\": " << uu << ", ";
    s << "\"R_Vy\": " << R_Vy << ", ";
    s << "\"R_Vu\": " << R_Vu << ", ";
    s << "\"R_umax\": " << R_umax << ", ";
    s << "\"alpha\": " << alpha << ", ";
    s << "\"gamma\": " << gamma << ", ";
    s << "\"beta\": " << beta << ", ";
    s << "\"symmetricMaxHistory\": " << (symmetricMaxHistory ? "true" : "false") << "}";
    return;
  }

  s << "MasonryShearMat tag: " << this->getTag() << endln;
  s << "  Ke: " << Ke << " initialVmax: " << initialVmax << " committedVmax: " << committedVmax << " Vmax: " << Vmax << " uu: " << uu << endln;
  s << "  R_Vy: " << R_Vy << " R_Vu: " << R_Vu << " R_umax: " << R_umax
    << " alpha: " << alpha << " gamma: " << gamma << " beta: " << beta << " symmetricMaxHistory: " << symmetricMaxHistory << endln;
  s << "  trial: strain=" << tState.strain << " stress=" << tState.stress << " tangent=" << tState.tangent << " branch=" << static_cast<int>(tState.branch) << " direction=" << tState.ldir << endln;
  s << "  committed: strain=" << cState.strain << " stress=" << cState.stress << " tangent=" << cState.tangent << " branch=" << static_cast<int>(cState.branch) << " direction=" << cState.ldir << endln;
  s << "  history: revStrain=" << cState.revStrain << " revStress=" << cState.revStress << " targetStrain=" << cState.tgtStrain << " targetStress=" << cState.tgtStress << " umaxPos=" << cState.umaxPos << " umaxNeg=" << cState.umaxNeg << endln;
}

// 创建只读诊断响应，供单元 recorder 或 OpenSeesPy 查询当前骨架与滞回状态。
Response *MasonryShearMat::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  if (argc > 0 && std::strcmp(argv[0], "diagnostics") == 0) {
    static const char *labels[] = {"Vmax", "committedVmax", "Vy", "uy", "branch", "strain", "stress", "tangent", "targetStrain", "targetStress", "reversalStrain", "reversalStress", "loadingDirection"};
    for (const char *label : labels) output.tag("ResponseType", label);
    return new MaterialResponse(this, 100, Vector(13));
  }
  return UniaxialMaterial::setResponse(argv, argc, output);
}

// 返回诊断向量，各分量顺序与 setResponse() 中的标签一致。
int MasonryShearMat::getResponse(int responseID, Information &information)
{
  if (responseID != 100) return UniaxialMaterial::getResponse(responseID, information);
  static Vector diagnostics(13);
  diagnostics(0) = Vmax;
  diagnostics(1) = committedVmax;
  diagnostics(2) = Vy;
  diagnostics(3) = uy;
  diagnostics(4) = static_cast<int>(tState.branch);
  diagnostics(5) = tState.strain;
  diagnostics(6) = tState.stress;
  diagnostics(7) = tState.tangent;
  diagnostics(8) = tState.tgtStrain;
  diagnostics(9) = tState.tgtStress;
  diagnostics(10) = tState.revStrain;
  diagnostics(11) = tState.revStress;
  diagnostics(12) = tState.ldir;
  return information.setVector(diagnostics);
}

// 根据当前试算最大剪切力更新所有依赖骨架峰值的派生参数。
void MasonryShearMat::updateDerivedParameters()
{
  Vy = R_Vy * Vmax;
  Vu = R_Vu * Vmax;
  umax = R_umax * uu;
  uy = Ke != 0.0 ? Vy / Ke : 0.0;
}

// 计算骨架曲线上的点
// u: 位移，可以为负值
void MasonryShearMat::backbone(double u, double &V, double &Et)
{
  if (std::abs(u) <= uy) {                       
    // 弹性段
    V = Ke * u;              
    Et = Ke;
  } else if (std::abs(u) <= umax) {              
    // 硬化段: 屈服点 -> 峰值点
    Et = (Vmax - Vy) / (umax - uy);
    double Vpos = Vy + Et * (std::abs(u) - uy);  
    V = u > 0 ? Vpos : -Vpos;  // 考虑 u 的符号
  } else {                        
    Et = (Vu - Vmax) / (uu - umax);   
    double Vpos = Vmax + Et * (std::abs(u) - umax);  
    Vpos = Vpos < 0 ? 0 : Vpos;  // 防止负值
    V = u > 0 ? Vpos : -Vpos;  // 考虑 u 的符号
  }
}

void MasonryShearMat::intersectionWithBackbone(double Strain, double Stress, double Tangent, double StrainDir, double &tgtStrain, double &tgtStress)
{
  // 目标是求出卸载/重加载直线
  //   sigma = Stress + Tangent * (u - Strain)
  // 与骨架曲线 backbone(u) 的交点。
  // 这里采用同向半轴上的二分搜索：只在当前应变的同号方向上寻找根，
  // 避免跨越零点引起的多个交点选择不稳定问题。

  const double eps = 1.0e-12;
  const double dir = (std::abs(StrainDir) > 1.0e-30) ? ((StrainDir > 0.0) ? 1.0 : -1.0) : ((Strain >= 0.0) ? 1.0 : -1.0);

  // 零应变时，交点即当前点，防止除零或无意义搜索。
  if (std::abs(Strain) <= eps) {
    tgtStrain = 0.0;
    tgtStress = 0.0;
    return;
  }

  // 确定搜索区间：限定在当前正/负半轴的历史峰值附近，
  // 因为卸载/重加载目标点通常落在同向骨架极值点上。
  double uLow = Strain;
  double uHigh = 0.0;

  if (dir > 0.0) {
    const double positiveHistory = cState.umaxPos > eps ? cState.umaxPos : -cState.umaxNeg;
    if (positiveHistory <= eps) {
      tgtStrain = Strain;
      tgtStress = Stress;
      return;
    }
    uHigh = std::max(positiveHistory, Strain + 1.0e-8);
  } else {
    const double negativeHistory = cState.umaxNeg < -eps ? cState.umaxNeg : -cState.umaxPos;
    if (std::abs(negativeHistory) <= eps) {
      tgtStrain = Strain;
      tgtStress = Stress;
      return;
    }
    uHigh = std::min(negativeHistory, Strain - 1.0e-8);
  }

  // 先检查当前点本身是否已在骨架曲线上。
  double Vb0 = 0.0;
  double tmpEt = 0.0;
  backbone(uLow, Vb0, tmpEt);
  const double fLow = Vb0 - (Stress + Tangent * (uLow - Strain));
  if (std::abs(fLow) <= 1.0e-10) {
    tgtStrain = uLow;
    tgtStress = Vb0;
    return;
  }

  double VbHigh = 0.0;
  backbone(uHigh, VbHigh, tmpEt);
  const double fHigh = VbHigh - (Stress + Tangent * (uHigh - Strain));

  // 如果当前区间不存在符号变化，退化为直接取同向骨架极值点。
  // 这是材料模型中常见的边界情形：交点落在极值附近或刚好与骨架相切。
  if (fLow * fHigh > 0.0) {
    tgtStrain = uHigh;
    tgtStress = VbHigh;
    return;
  }

  double uLeft = uLow;
  double uRight = uHigh;
  double VbLeft = Vb0;
  double fLeft = fLow;
  double uMid = 0.0;
  double VbMid = 0.0;
  double fMid = 0.0;

  // 二分法查找交点：在 [uLow, uHigh] 上求解骨架曲线与卸载/重加载射线的根。
  for (int iter = 0; iter < 80; ++iter) {
    uMid = 0.5 * (uLeft + uRight);
    backbone(uMid, VbMid, tmpEt);
    fMid = VbMid - (Stress + Tangent * (uMid - Strain));

    if (std::abs(fMid) <= 1.0e-10 || std::abs(uRight - uLeft) <= 1.0e-12) {
      tgtStrain = uMid;
      tgtStress = VbMid;
      return;
    }

    if (fLeft * fMid <= 0.0) {
      uRight = uMid;
    } else {
      uLeft = uMid;
      fLeft = fMid;
      VbLeft = VbMid;
    }
  }

  // 若迭代未完全收敛，返回最后一次估计值，保证调用方仍然有有效目标点。
  tgtStrain = uMid;
  tgtStress = VbMid;
}

// 返回状态清零的 MasonryShearMat::State, 仅切线刚度为 Ke
MasonryShearMat::State
MasonryShearMat::initialState() const
{
    State s = {};            // 全部清零
    s.tangent = Ke;          // 初始切线刚度
    return s;
}

// 更新试状态是否变化
bool MasonryShearMat::checkTransitions() {

  if (fabs(tState.strain - cState.strain) < 1.0e-14) return false;

  // 从 ELASTIC 状态变化
  if (cState.branch == ELASTIC) {
    if (tState.strain > uy) {
      // 应变超过屈服位移, 进入 HARDENING_1 状态
      tState.branch = HARDENING_1;
      return true;
    } else if (tState.strain < -uy) {
      // 应变超过负向屈服位移, 进入 HARDENING_1 状态
      tState.branch = HARDENING_1;
      return true;
    } else {
      // 保持在 ELASTIC 状态
      return false;
    }
  } 

  // 从 HARDENING_1 状态变化
  if (cState.branch == HARDENING_1) {
    if (tState.strain > umax || tState.strain < -umax) {
      // 应变超过最大位移, 进入 HARDENING_2 状态
      tState.branch = HARDENING_2;
      return true;
    } else if (tState.ldir * cState.strain < 0.0) {
      // 卸载
      tState.branch = UNLOADING_1;
      return true;
    } else {
      // 保持在 HARDENING_1 状态
      return false;
    }
  }

  // 从 HARDENING_2 状态变化
  if (cState.branch == HARDENING_2) {
    if (tState.ldir * cState.strain < 0.0) {
      // 卸载
      tState.branch = UNLOADING_1;
      return true;
    } else {
      // 保持在 HARDENING_2 状态
      return false; 
    }
  }

  // 从 UNLOADING_1 状态变化
  if (cState.branch == UNLOADING_1) {
    if (tState.ldir * cState.stress > 0.0) {
      // 反向加载
      tState.branch = RELOADING_FROM_UNLOADING_1;
      return true;
    } else {
      // 正向加载
      // 计算假设还在 UNLOADING_1 状态下的应力
      double tmpStress = cState.tangent * (tState.strain - cState.strain) + cState.stress;
      if (fabs(tmpStress) < gamma * fabs(cState.revStress)) {
        // 进入 UNLOADING_2 状态
        tState.branch = UNLOADING_2;
        return true;
      } else {
        // 保持在 UNLOADING_1 状态
        tState.branch = UNLOADING_1;
        return false;
      }
    } 
  }

  // 从 UNLOADING_2 状态变化
  if (cState.branch == UNLOADING_2) {
    if (tState.ldir * cState.ldir < 0.0) {
      // 反向加载
        if ( (tState.ldir > 0 && cState.stress < gamma * cState.tgtStress) || (tState.ldir < 0 && cState.stress > gamma * cState.tgtStress) ) {
        // 反向加载到 UNLOADING_1 状态
        tState.branch = UNLOADING_1;
      } else {
        // 反向加载到 UNLOADING_2 状态
        tState.branch = UNLOADING_2;
      }
      return true;
    } else if ((tState.strain - tState.tgtStrain) * tState.ldir > 0.0) {
      // 进入 HARDENING 状态
      if (std::abs(tState.strain) > umax) {
        tState.branch = HARDENING_2;
      } else {
        tState.branch = HARDENING_1;
      }
      return true;
    } else {
      // 保持在 UNLOADING_2 状态
      return false;
    }
  }

  // 从 RELOADING_FROM_UNLOADING_1 状态变化
  if (cState.branch == RELOADING_FROM_UNLOADING_1) {
    if (tState.ldir * cState.stress < 0.0) {
      // 局部加载方向反转时优先进入卸载，不能因跨过旧目标点而误回骨架。
      tState.branch = UNLOADING_1;
      return true;
    } else if ((tState.strain - cState.tgtStrain) * tState.ldir > 0.0) {
      // 进入 HARDENING 状态
      if (std::abs(tState.strain) > umax) {
        tState.branch = HARDENING_2;
      } else {
        tState.branch = HARDENING_1;
      }
      return true;
    } else {
      // 还是处于 RELOADING_FROM_UNLOADING_1 状态
      return false;
    }
  }

  return false;
}

void MasonryShearMat::ruleFromElastic() {
  if (tState.branch == ELASTIC) {
    // 更新应力
    tState.stress = Ke * tState.strain;
    tState.tangent = Ke;
  } else if (tState.branch == HARDENING_1) {
    // 更新应力和切线刚度
    backbone(tState.strain, tState.stress, tState.tangent);
  } 
}

void MasonryShearMat::ruleFromHardening1() {
  if (tState.branch == HARDENING_1) {
    // 更新应力和切线刚度
    backbone(tState.strain, tState.stress, tState.tangent);
  } else if (tState.branch == UNLOADING_1) {
    // 更新应力和切线刚度
    tState.tangent = computeUnload1Tangent(cState);
    tState.stress = tState.tangent * (tState.strain - cState.strain) + cState.stress;
    // 从骨架反转时开启新的耗能回路。
    tState.cycleStartStrain = cState.strain;
    tState.cycleStartStress = cState.stress;
    tState.cycleStartDir = cState.strain >= 0.0 ? 1.0 : -1.0;
    tState.cycleWork = 0.0;
    // 更新反转点历史变量
    tState.revStrain = cState.strain;
    tState.revStress = cState.stress;
  } else if (tState.branch == HARDENING_2) {
    // 更新应力和切线刚度
    backbone(tState.strain, tState.stress, tState.tangent);
  }
}

void MasonryShearMat::ruleFromHardening2() {
  if (tState.branch == HARDENING_2) {
    // 更新应力和切线刚度
    backbone(tState.strain, tState.stress, tState.tangent);
  } else if (tState.branch == UNLOADING_1) {
    // 更新应力和切线刚度
    tState.tangent = computeUnload1Tangent(cState);
    tState.stress = tState.tangent * (tState.strain - cState.strain) + cState.stress;
    // 从骨架反转时开启新的耗能回路。
    tState.cycleStartStrain = cState.strain;
    tState.cycleStartStress = cState.stress;
    tState.cycleStartDir = cState.strain >= 0.0 ? 1.0 : -1.0;
    tState.cycleWork = 0.0;
    // 更新反转点历史变量
    tState.revStrain = cState.strain;
    tState.revStress = cState.stress;
  }
}

void MasonryShearMat::ruleFromUnloading1() {
  if (tState.branch == UNLOADING_1) {
    // 更新应力
    tState.stress = tState.tangent * (tState.strain - cState.strain) + cState.stress;
  } else if (tState.branch == RELOADING_FROM_UNLOADING_1) {
    // 更新应力
    tState.stress = tState.tangent * (tState.strain - cState.strain) + cState.stress;
    // 更新目标点
    const double strainDir = (tState.strain - cState.strain) >= 0.0 ? 1.0 : -1.0;
    intersectionWithBackbone(cState.strain, cState.stress, cState.tangent, strainDir, tState.tgtStrain, tState.tgtStress);
  } else if (tState.branch == UNLOADING_2) {
    // 以未退化历史最大位移闭合耗能回路，再设置强度退化后的目标点。
    setUnloading2Target(tState, cState);
    // 更新应力和切线刚度
    tState.tangent = (tState.tgtStress - cState.stress) / (tState.tgtStrain - cState.strain);
    tState.stress = cState.stress + tState.tangent * (tState.strain - cState.strain);
  }
}

void MasonryShearMat::ruleFromUnloading2() {
  if (tState.branch == UNLOADING_2) {
    // 有两种可能，一种继续加载到目标点，另一种是反向加载
    const double du = tState.strain - cState.strain;
    if (std::abs(du) <= 1.0e-14) {
      // 零增量只保持当前状态，不能被识别为反向加载
      tState.stress = cState.stress;
    } else if (du * cState.tgtStrain > 0.0) {
      // 继续加载到目标点
      // 更新应力
      tState.stress = cState.stress + tState.tangent * (tState.strain - cState.strain);
    } else {
      // 反向加载
      // 更新应力和切线刚度
      tState.tangent = computeUnload1Tangent(cState);
      tState.stress = cState.stress + tState.tangent * (tState.strain - cState.strain);
      // 更新反转点
      tState.revStrain = cState.strain;
      tState.revStress = cState.stress;
      // 以未退化历史最大位移闭合耗能回路，再设置强度退化后的目标点。
      setUnloading2Target(tState, cState);
    }
  } else if (tState.branch == HARDENING_1 || tState.branch == HARDENING_2) {
    // 更新应力和切线刚度
    backbone(tState.strain, tState.stress, tState.tangent);
  } else if (tState.branch == UNLOADING_1) {
    // 反向加载
    // 更新应力和切线刚度
    tState.tangent = computeUnload1Tangent(cState);
    tState.stress = cState.stress + tState.tangent * (tState.strain - cState.strain);
    // 更新反转点
    tState.revStrain = cState.strain;
    tState.revStress = cState.stress;
  }
}

void MasonryShearMat::ruleFromReloadingFromUnloading1() {
  if (tState.branch == RELOADING_FROM_UNLOADING_1) {
    // 更新应力
    tState.stress = cState.tangent * (tState.strain - cState.strain) + cState.stress;
  } else if (tState.branch == HARDENING_1 || tState.branch == HARDENING_2) {
    // 更新应力和切线刚度
    backbone(tState.strain, tState.stress, tState.tangent);
  } else if (tState.branch == UNLOADING_1) {
    // 以当前反转点重新计算第一卸载段刚度和应力。
    tState.tangent = computeUnload1Tangent(cState);
    tState.stress = tState.tangent * (tState.strain - cState.strain) + cState.stress;
  }
}

double MasonryShearMat::computeUnload1Tangent(const State &state) {
  // 根据当前反转一侧的历史最大位移计算卸载刚度。
  const double maximumStrain = state.strain >= 0.0 ? state.umaxPos : -state.umaxNeg;
  if (maximumStrain <= uy)
    return Ke;

  const double Ck = (alpha - 1.0) / (uu / uy - 1.0);
  return Ke * (1.0 + Ck * (maximumStrain / uy - 1.0));
}

double MasonryShearMat::backboneWork(double strain) const {
  // 骨架曲线关于原点反对称，因此从零点到正、负同幅值位移的功相同。
  const double absoluteStrain = std::fabs(strain);
  if (absoluteStrain <= uy)
    return 0.5 * Ke * absoluteStrain * absoluteStrain;

  const double workAtYield = 0.5 * Vy * uy;
  const double hardening1Tangent = (Vmax - Vy) / (umax - uy);
  if (absoluteStrain <= umax) {
    const double du1 = absoluteStrain - uy;
    return workAtYield + Vy * du1 + 0.5 * hardening1Tangent * du1 * du1;
  }

  const double hardening1Length = umax - uy;
  const double workAtPeak = workAtYield + Vy * hardening1Length + 0.5 * hardening1Tangent * hardening1Length * hardening1Length;
  const double hardening2Tangent = (Vu - Vmax) / (uu - umax);
  double hardening2Length = absoluteStrain - umax;
  if (hardening2Tangent < 0.0)
    hardening2Length = std::min(hardening2Length, -Vmax / hardening2Tangent);
  return workAtPeak + Vmax * hardening2Length + 0.5 * hardening2Tangent * hardening2Length * hardening2Length;
}

void MasonryShearMat::setUnloading2Target(State &state, const State &origin) {
  // 优先使用目标方向自己的历史峰值；该方向尚未加载时，仅对重加载目标点镜像另一方向历史。
  double originalTargetStrain = state.ldir > 0.0 ? origin.umaxPos : origin.umaxNeg;
  if (state.ldir > 0.0 && originalTargetStrain <= 0.0)
    originalTargetStrain = -origin.umaxNeg;
  else if (state.ldir < 0.0 && originalTargetStrain >= 0.0)
    originalTargetStrain = -origin.umaxPos;
  double originalTargetStress = 0.0;
  double originalTargetTangent = 0.0;
  backbone(originalTargetStrain, originalTargetStress, originalTargetTangent);

  // 用当前路径、通向原目标的第二卸载段以及返回回路起点的骨架段闭合回路。
  double cycleEnergy = 0.0;
  if (origin.cycleStartDir != 0.0) {
    const double projectedWork = 0.5 * (origin.stress + originalTargetStress) * (originalTargetStrain - origin.strain);
    const double closingBackboneWork = backboneWork(origin.cycleStartStrain) - backboneWork(originalTargetStrain);
    cycleEnergy = std::fabs(origin.cycleWork + projectedWork + closingBackboneWork);
  }

  // 峰值以前向外移动目标会提高强度，因此强度退化只在原目标达到峰值位移后启用。
  const double degradationIncrement = std::fabs(originalTargetStrain) >= umax ? beta * cycleEnergy / Vmax : 0.0;
  state.tgtStrain = originalTargetStrain + (originalTargetStrain >= 0.0 ? degradationIncrement : -degradationIncrement);
  double targetTangent = 0.0;
  backbone(state.tgtStrain, state.tgtStress, targetTangent);
}

void MasonryShearMat::accumulateCycleWork() {
  // 未从骨架开启耗能回路时不累计单调加载功。
  if (tState.cycleStartDir == 0.0)
    return;

  const double strainIncrement = tState.strain - cState.strain;
  if (std::fabs(strainIncrement) > 1.0e-14)
    tState.cycleWork += 0.5 * (cState.stress + tState.stress) * strainIncrement;

  // 回到骨架后结束当前耗能回路并清除临时累计量。
  if (tState.branch == HARDENING_1 || tState.branch == HARDENING_2 || tState.branch == ELASTIC) {
    tState.cycleStartStrain = 0.0;
    tState.cycleStartStress = 0.0;
    tState.cycleStartDir = 0.0;
    tState.cycleWork = 0.0;
  }
}
