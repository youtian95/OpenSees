// 模块功能：定义三维砌体宏单元 MasonryMacro3D 的 OpenSees 接口与数据结构。
// 使用流程：通过 element masonryMacro 命令创建单元，框架负责节点、材料、坐标变换、状态和并行通讯。
// 输入输出：输入两个六自由度节点、六类单轴材料和坐标变换；输出 12x12 刚度与 12 维恢复力。

#ifndef MasonryMacro3D_h
#define MasonryMacro3D_h

#include <Element.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>

class Channel;
class CrdTransf;
class Domain;
class FEM_ObjectBroker;
class Information;
class Node;
class OPS_Stream;
class Renderer;
class Response;
class UniaxialMaterial;

// 三维两节点砌体宏单元：每个节点具有 ux、uy、uz、rx、ry、rz 六个自由度。
class MasonryMacro3D : public Element
{
public:
    // 创建供分析使用的三维砌体宏单元。
    // tag: 单元标签。
    // nodeI: I 端节点标签。
    // nodeJ: J 端节点标签。
    // bendingYMaterial: 绕局部 y 轴弯曲材料原型，I、J 两端各复制一份。
    // bendingZMaterial: 绕局部 z 轴弯曲材料原型，I、J 两端各复制一份。
    // shearYMaterial: 局部 y 方向剪切材料原型。
    // shearZMaterial: 局部 z 方向剪切材料原型。
    // torsionMaterial: 绕局部 x 轴扭转材料原型。
    // axialMaterial: 局部 x 方向轴向材料原型。
    // coordinateTransformation: 三维坐标变换原型，局部轴方向由 geomTransf 定义。
    MasonryMacro3D(int tag, int nodeI, int nodeJ, UniaxialMaterial &bendingYMaterial, UniaxialMaterial &bendingZMaterial, UniaxialMaterial &shearYMaterial, UniaxialMaterial &shearZMaterial, UniaxialMaterial &torsionMaterial, UniaxialMaterial &axialMaterial, CrdTransf &coordinateTransformation);

    // 创建供 ObjectBroker 接收数据使用的空对象。
    MasonryMacro3D();

    // 释放单元拥有的材料副本和坐标变换副本。
    ~MasonryMacro3D();

    const char *getClassType(void) const { return "MasonryMacro3D"; }

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);
    int getNumDOF(void);

    // theDomain: 单元所属的 Domain；传入空指针表示解除关联。
    void setDomain(Domain *theDomain);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    int update(void);

    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Vector &getResistingForce(void);

    // commitTag: 当前提交标签。
    // theChannel: 接收序列化数据的通道。
    int sendSelf(int commitTag, Channel &theChannel);

    // commitTag: 当前提交标签。
    // theChannel: 提供序列化数据的通道。
    // theBroker: 根据 classTag 创建材料和坐标变换的对象工厂。
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    // theViewer: 图形渲染器。
    // displayMode: 显示模式。
    // fact: 位移放大系数。
    // modes: 可选显示模式名称数组。
    // numModes: modes 数组长度。
    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes = 0, int numModes = 0);

    // s: 输出流。
    // flag: 输出格式标志。
    void Print(OPS_Stream &s, int flag = 0);

    // argv: 响应请求参数。
    // argc: 响应请求参数数量。
    // output: 响应描述输出流。
    Response *setResponse(const char **argv, int argc, OPS_Stream &output);

    // responseID: setResponse 分配的响应编号。
    // information: 承载响应结果的对象。
    int getResponse(int responseID, Information &information);

private:
    ID connectedExternalNodes;
    Node *theNodes[2];
    CrdTransf *theCoordTransf;

    UniaxialMaterial *bendingYMaterials[2];
    UniaxialMaterial *bendingZMaterials[2];
    UniaxialMaterial *shearYMaterial;
    UniaxialMaterial *shearZMaterial;
    UniaxialMaterial *torsionMaterial;
    UniaxialMaterial *axialMaterial;

    // 以下对象作为 OpenSees const 引用返回值的持久缓存。
    Matrix tangentStiffness;
    Matrix initialStiffness;
    Vector resistingForce;
};

#endif
