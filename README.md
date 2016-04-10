# GraduateProject
Leap Motion + BulletPhysics + OpenGL

我在Bullet這個開源的物理引擎內增加了一種新的形變，Bullet本身沒有撕裂的效果，只有在柔性物體上打洞的效果，利用Bullet內柔性物體的
組成特性，去製造物體破碎的效果。

選用的柔體是由三角形組成，我做的事情主要就是根據使用者剪到的位置不斷去新增新的節點，不斷的打破原本的結構
，在不斷的重新連成三角形。
