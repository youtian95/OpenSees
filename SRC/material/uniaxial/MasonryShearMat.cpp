// ============================================================================
// 文件: MasonryShearMat.cpp
// 功能: MasonryShearMat 砌体剪切单轴材料实现 + OPS 工厂函数
// 说明: 样板代码(工厂函数/状态管理/并行通讯)已完成,
//       核心本构规律在 setTrialStrain 中以 TODO 标出, 由开发者手动实现
// 用法: 编译后在脚本中使用:
//       Tcl:        uniaxialMaterial MasonryShearMat matTag? Ke? Vmax? uu? <可选参数>
//       OpenSeesPy: ops.uniaxialMaterial('MasonryShearMat', matTag, Ke, Vmax, uu, ...)
// 输入: matTag(整数); Ke Vmax uu(必选实数); R_Vy~beta(可选, 尾部省略用默认值)
// 输出: 创建 MasonryShearMat 材料对象并注册到模型构建器
// ============================================================================

#include "MasonryShearMat.h"

#include <OPS_Globals.h>
#include <elementAPI.h>
#include <classTags.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Vector.h>
#include <cmath>

// 已创建的材料实例数(用于只在首次创建时打印署名横幅)
static int numMasonryShearMat = 0;

// 工厂函数: uniaxialMaterial 命令的入口
// 命令语法(< > 内为可选参数, 尾部省略时使用构造函数的默认参数):
//   uniaxialMaterial MasonryShearMat matTag? Ke? Vmax? uu? <R_Vy?> <R_Vu?> <R_umax?> <alpha?> <gamma?> <beta?>
void *OPS_MasonryShearMat(void)
{
  // 首次创建时打印一次署名信息
  if (numMasonryShearMat == 0) {
    opserr << "MasonryShearMat uniaxial material - Written by youtian95\n";
    numMasonryShearMat = 1;
  }

  // 至少需要 matTag + Ke + Vmax + uu 四个参数
  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "WARNING invalid args, want: uniaxialMaterial MasonryShearMat tag? Ke? Vmax? uu? <R_Vy?> <R_Vu?> <R_umax?> <alpha?> <gamma?> <beta?>\n";
    return 0;
  }

  // 1) 读取整数参数: matTag
  int iData[1];
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial MasonryShearMat tag" << endln;
    return 0;
  }

  // 2) 读取实数参数: 只读用户实际给出的个数(封顶9个),
  //    再按实际个数调用构造函数, 尾部未给出的参数由头文件中的默认参数补齐
  double dData[9];
  int numArgs = OPS_GetNumRemainingInputArgs();
  numData = (numArgs > 9) ? 9 : numArgs;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid MasonryShearMat parameters\n";
    return 0;
  }

  // 3) 按实际读到的参数个数创建材料, 缺省参数自动生效
  UniaxialMaterial *theMaterial = 0;
  switch (numData) {
    case 3:  theMaterial = new MasonryShearMat(iData[0], dData[0], dData[1], dData[2]); break;
    case 4:  theMaterial = new MasonryShearMat(iData[0], dData[0], dData[1], dData[2], dData[3]); break;
    case 5:  theMaterial = new MasonryShearMat(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4]); break;
    case 6:  theMaterial = new MasonryShearMat(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]); break;
    case 7:  theMaterial = new MasonryShearMat(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6]); break;
    case 8:  theMaterial = new MasonryShearMat(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7]); break;
    default: theMaterial = new MasonryShearMat(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], dData[8]); break;
  }

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type MasonryShearMat\n";
    return 0;
  }
  return theMaterial;
}

// 全参构造函数
// tag: 材料标签; 其余参数含义见头文件
// 基类必须用 (tag, MAT_TAG_MasonryShearMat) 初始化
MasonryShearMat::MasonryShearMat(int tag, double _Ke, double _Vmax, double _uu, double _R_Vy, double _R_Vu, double _R_umax, double _alpha, double _gamma, double _beta)
  : UniaxialMaterial(tag, MAT_TAG_MasonryShearMat),
    Ke(_Ke), Vmax(_Vmax), uu(_uu),
    R_Vy(_R_Vy), R_Vu(_R_Vu), R_umax(_R_umax),
    alpha(_alpha), gamma(_gamma), beta(_beta)
{
  // TODO: 如有其他派生参数/历史变量, 在此初始化
  // 派生参数计算: 屈服力 Vy, 极限位移对应剪力 Vu, 最大剪切力对应位移 umax
  Vy = R_Vy * Vmax;
  Vu = R_Vu * Vmax;
  umax = R_umax * uu;
  uy = Vy / Ke;
  // 历史变量初始化为零
  tState = initialState();
  cState = initialState();
}

// 默认构造函数: 供 FEM_ObjectBroker 在并行/数据库恢复时使用
MasonryShearMat::MasonryShearMat()
  : UniaxialMaterial(0, MAT_TAG_MasonryShearMat),
    Ke(0.0), Vmax(0.0), uu(0.0),
    R_Vy(0.7), R_Vu(0.8), R_umax(0.5),
    alpha(0.8), gamma(0.6), beta(0.06),
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

  return 0;
}

// 以下四个 get 函数返回当前试状态
double MasonryShearMat::getStrain(void)         { return tState.strain; }
double MasonryShearMat::getStress(void)         { return tState.stress; }
double MasonryShearMat::getTangent(void)        { return tState.tangent; }
double MasonryShearMat::getInitialTangent(void) { return Ke; }

// 提交状态: 试状态转为已提交状态
int
MasonryShearMat::commitState(void)
{
  cState = tState;
  return 0;
}

// 回滚到上次提交状态(迭代不收敛时调用)
int
MasonryShearMat::revertToLastCommit(void)
{
  tState = cState;
  return 0;
}

// 回滚到初始状态(重新开始分析时调用)
int
MasonryShearMat::revertToStart(void)
{
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

// 打包发送(并行/数据库): 材料tag + 9个材料参数 + 完整已提交状态
// commitTag: 提交标签; theChannel: 通讯通道
int
MasonryShearMat::sendSelf(int commitTag, Channel &theChannel)
{
  // 打包布局: data(0)=材料tag; data(1..9)=材料参数; data(10..21)=已提交状态全字段
  static Vector data(22);
  data(0) = this->getTag();
  // 材料参数
  data(1) = Ke;     data(2) = Vmax;   data(3) = uu;
  data(4) = R_Vy;   data(5) = R_Vu;   data(6) = R_umax;
  data(7) = alpha;  data(8) = gamma;  data(9) = beta;
  // 已提交状态: 基本量 + 当前分支 + 全部历史变量
  data(10) = cState.strain;
  data(11) = cState.strainRate;
  data(12) = cState.stress;
  data(13) = cState.tangent;
  data(14) = (int)cState.branch;
  data(15) = cState.ldir;
  data(16) = cState.umaxPos;
  data(17) = cState.umaxNeg;
  data(18) = cState.revStrain;
  data(19) = cState.revStress;
  data(20) = cState.tgtStrain;
  data(21) = cState.tgtStress;

  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0)
    opserr << "MasonryShearMat::sendSelf() - failed to send data\n";
  return res;
}

// 接收恢复(并行/数据库): 打包顺序必须与 sendSelf 严格对应
int
MasonryShearMat::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static Vector data(22);
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "MasonryShearMat::recvSelf() - failed to receive data\n";
    return res;
  }
  this->setTag((int)data(0));
  // 材料参数
  Ke = data(1);     Vmax = data(2);   uu = data(3);
  R_Vy = data(4);   R_Vu = data(5);   R_umax = data(6);
  alpha = data(7);  gamma = data(8);  beta = data(9);
  // 重算派生参数(派生量不传输, 必须由参数重建)
  Vy = R_Vy * Vmax;
  Vu = R_Vu * Vmax;
  umax = R_umax * uu;
  uy = (Ke > 0.0) ? Vy / Ke : 0.0;
  // 已提交状态: 基本量 + 当前分支 + 全部历史变量
  cState.strain     = data(10);
  cState.strainRate = data(11);
  cState.stress     = data(12);
  cState.tangent    = data(13);
  cState.branch     = (Branch)(int)data(14);
  cState.ldir       = data(15);
  cState.umaxPos    = data(16);
  cState.umaxNeg    = data(17);
  cState.revStrain  = data(18);
  cState.revStress  = data(19);
  cState.tgtStrain  = data(20);
  cState.tgtStress  = data(21);
  // 用已提交状态初始化试状态
  tState = cState;
  return 0;
}

// 打印材料信息
// s: 输出流; flag: 输出详细程度标志
void
MasonryShearMat::Print(OPS_Stream &s, int flag)
{
  s << "MasonryShearMat tag: " << this->getTag() << endln;
  s << "  Ke: " << Ke << " Vmax: " << Vmax << " uu: " << uu << endln;
  s << "  R_Vy: " << R_Vy << " R_Vu: " << R_Vu << " R_umax: " << R_umax
    << " alpha: " << alpha << " gamma: " << gamma << " beta: " << beta << endln;
  // TODO: 按需打印更多状态信息
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
    if (cState.umaxPos <= eps) {
      tgtStrain = Strain;
      tgtStress = Stress;
      return;
    }
    uHigh = std::max(cState.umaxPos, Strain + 1.0e-8);
  } else {
    if (std::abs(cState.umaxNeg) <= eps) {
      tgtStrain = Strain;
      tgtStress = Stress;
      return;
    }
    uHigh = std::min(cState.umaxNeg, Strain - 1.0e-8);
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
    if ((tState.strain - cState.tgtStrain) * tState.ldir > 0.0) {
      // 进入 HARDENING 状态
      if (std::abs(tState.strain) > umax) {
        tState.branch = HARDENING_2;
      } else {
        tState.branch = HARDENING_1;
      }
      return true;
    } else if ( tState.ldir * cState.stress < 0.0) {
      // 进入 UNLOADING_1 状态
      tState.branch = UNLOADING_1;
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
    // 更新历史变量
    tState.umaxPos = std::max(tState.umaxPos, tState.strain);
    tState.umaxNeg = std::min(tState.umaxNeg, tState.strain);
  } else if (tState.branch == HARDENING_1) {
    // 更新应力和切线刚度
    backbone(tState.strain, tState.stress, tState.tangent);
    // 更新历史变量
    tState.umaxPos = std::max(std::max(tState.umaxPos, tState.strain), uy);
    tState.umaxNeg = std::min(std::min(tState.umaxNeg, tState.strain), -uy);
  } 
}

void MasonryShearMat::ruleFromHardening1() {
  if (tState.branch == HARDENING_1) {
    // 更新应力和切线刚度
    backbone(tState.strain, tState.stress, tState.tangent);
    // 更新历史变量
    tState.umaxPos = std::max(tState.umaxPos, tState.strain);
    tState.umaxNeg = std::min(tState.umaxNeg, tState.strain);
  } else if (tState.branch == UNLOADING_1) {
    // 更新应力和切线刚度
    tState.tangent = computeUnload1Tangent(cState);
    tState.stress = tState.tangent * (tState.strain - cState.strain) + cState.stress;
    // 更新反转点历史变量
    tState.revStrain = cState.strain;
    tState.revStress = cState.stress;
  } else if (tState.branch == HARDENING_2) {
    // 更新应力和切线刚度
    backbone(tState.strain, tState.stress, tState.tangent);
    // 更新历史变量
    tState.umaxPos = std::max(tState.umaxPos, tState.strain);
    tState.umaxNeg = std::min(tState.umaxNeg, tState.strain);
  }
}

void MasonryShearMat::ruleFromHardening2() {
  if (tState.branch == HARDENING_2) {
    // 更新应力和切线刚度
    backbone(tState.strain, tState.stress, tState.tangent);
    // 更新历史变量
    tState.umaxPos = std::max(tState.umaxPos, tState.strain);
    tState.umaxNeg = std::min(tState.umaxNeg, tState.strain);
  } else if (tState.branch == UNLOADING_1) {
    // 更新应力和切线刚度
    tState.tangent = computeUnload1Tangent(cState);
    tState.stress = tState.tangent * (tState.strain - cState.strain) + cState.stress;
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
    // 更新目标点
    double targetTangent;
    if ((tState.strain - cState.strain) < 0.0) {
      tState.tgtStrain = cState.umaxNeg;
      backbone(tState.tgtStrain, tState.tgtStress, targetTangent);
    } else {
      tState.tgtStrain = cState.umaxPos;
      backbone(tState.tgtStrain, tState.tgtStress, targetTangent);
    }
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
      // 更新目标点
      if ((tState.strain - cState.strain) < 0.0) {
        tState.tgtStrain = cState.umaxNeg;
        double targetTangent;
        backbone(tState.tgtStrain, tState.tgtStress, targetTangent);
      } else {
        tState.tgtStrain = cState.umaxPos;
        double targetTangent;
        backbone(tState.tgtStrain, tState.tgtStress, targetTangent);
      }
    }
  } else if (tState.branch == HARDENING_1 || tState.branch == HARDENING_2) {
    // 更新应力和切线刚度
    backbone(tState.strain, tState.stress, tState.tangent);
    // 更新历史变量
    tState.umaxPos = std::max(tState.umaxPos, tState.strain);
    tState.umaxNeg = std::min(tState.umaxNeg, tState.strain);
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
    // 更新历史变量
    tState.umaxPos = std::max(tState.umaxPos, tState.strain);
    tState.umaxNeg = std::min(tState.umaxNeg, tState.strain);
  } else if (tState.branch == UNLOADING_1) {
    // 以当前反转点重新计算第一卸载段刚度和应力。
    tState.tangent = computeUnload1Tangent(cState);
    tState.stress = tState.tangent * (tState.strain - cState.strain) + cState.stress;
  }
}

double MasonryShearMat::computeUnload1Tangent(const State &state) {
  // 在屈服点取Ke，在极限位移点取alpha*Ke，两个骨架硬化阶段使用同一条连续关系。
  const double Ck = (alpha - 1.0) / (uu / uy - 1.0);
  return Ke * (1.0 + Ck * (std::fabs(state.strain) / uy - 1.0));
}
