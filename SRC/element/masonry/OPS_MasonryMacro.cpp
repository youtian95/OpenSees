// 模块功能：为 masonryMacro 命令提供二维/三维统一入口。
// 使用流程：OpenSees 解析 element masonryMacro 后调用本函数，并依据当前 NDM 分派到对应实现。
// 输入输出：输入当前模型维数和剩余命令参数；输出新建的 MasonryMacro2D 或 MasonryMacro3D 指针。

#include <OPS_Globals.h>
#include <elementAPI.h>

void *OPS_MasonryMacro2D(void);
void *OPS_MasonryMacro3D(void);

void *OPS_MasonryMacro(void)
{
    int modelDimension = OPS_GetNDM();
    if (modelDimension == 2) {
        return OPS_MasonryMacro2D();
    }
    if (modelDimension == 3) {
        return OPS_MasonryMacro3D();
    }

    opserr << "WARNING masonryMacro supports only 2D and 3D models" << endln;
    return 0;
}
