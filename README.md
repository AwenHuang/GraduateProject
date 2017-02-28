# GraduateProject
Leap Motion + BulletPhysics + OpenGL
Visual Studio 2010/2012

DEMO Video:
https://youtu.be/ymW3nfn7wBQ

用Leap Motion偵測手勢，改寫Lepa Motion的API新增一個剪刀手的手勢的判斷，需要注意的是Leap Motion的Z軸剛好跟Bullet相反

OpenGL建構場景還有紅藍效果

我們在Bullet這個開源的物理引擎內增加了一種新的形變，Bullet本身沒有撕裂的效果，只有在柔性物體上打洞的效果，利用Bullet內柔性物體的
組成特性，去製造物體破碎的效果。

選用的柔體是由三角形組成，我做的事情主要就是根據使用者剪到的位置不斷去新增新的節點，不斷的打破原本的結構
，在不斷的重新連成三角形。

形變的主體code主要寫在
src -> BulletSoftBody -> btSoftBody.cpp, 1558行開始


