This is my attempt on COMP5821M Geometric Processing from University of Leeds

Video demonstration:  
This code can run on University machine but have not tested for other platforms yet.  
  
Video demonstration: https://www.youtube.com/watch?v=eF6QlOhJGGo


To compile on feng-linux / feng-gps:

module add qt/5.13.0
qmake -project QT+=opengl
qmake
make

To compile on OSX:
Use Homebrew to install qt

qmake -project QT+=opengl
qmake
make

To compile on Windows:
Unfortunately, the official OpenGL on Windows was locked at GL 1.1.  Many many hacks exist, and they all disagree.
Just to make it worse, the Qt response to this is clumsy.  Net result: there is no easy way to get this compiling on Windows.
I will aim to update these instructions at a later date.

To run on feng-linux / feng-gps:

./LoopSubdivisionRelease ../path_to/model.diredgenormal

To run on OSX:
./LoopSubdivisionRelease.app/Contents/MacOS/FakeGLRenderWindowRelease  ../path_to/model.diredgenormal

