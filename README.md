360° video player
======

The 360° video player supports diverse projection schemes, such as equi-rectangular, adjusted equal-area, standard cubemap, and equi-angular cubemap. 
We implement the 360° video player in C++ language, which is used to display test video sequences to HMD.
The player is also connected with OSVR server to store the subjects’ head movement data aligned with each frame.

Prerequisites
---------------

- libav library from ffmpeg. (I currently use FFmpeg 3.3)

- [Glew-2.1.0-win32](https://sourceforge.net/projects/glew/files/glew/2.1.0/)

- [boost_1_57_0](https://www.boost.org/users/history/version_1_57_0.html) with Python 3.4

- [libzmq 4.2.2](https://github.com/zeromq/libzmq/releases)

- [SDL2](https://www.libsdl.org/download-2.0.php)

- [SOIL](https://github.com/kbranigan/Simple-OpenGL-Image-Library)

- [OSVR-Core](https://github.com/OSVR/OSVR-Core) and [OSVR-Rendermanager](https://github.com/sensics/OSVR-RenderManager)
libraries from the OSVR opensource project. 


Compilation steps
-----------------

mkdir build > cd build > cmake -G "Visual Studio 15 2017" ..

Tested platform:
----------------

We used the Oculus DK2 HMD, which is set up on a PC with Intel i7-3770 CPU, NVIDIA GTX
1060 GPU, 16 GB RAM, and Windows 10 OS.

Acknowledgments
----------------

Thanks for the [contributions](https://github.com/xmar/360Degree_Head_Movement_Dataset) proposed by Corbillon et al., which is enhanced by us for additional projection schemes.

Contacts:
---------
E-mail: shunhuai.yao@gapp.nthu.edu.tw
