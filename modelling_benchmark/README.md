A introduction Video is accessible from https://youtu.be/hJ_x5sJbg1w.

Do not update this folder as it is benchmarked by the youtube video

update on 240808

it is found that msys64 is quite straightforward to install gfortran and git in windows enviroment, as it uses pacman to managet packages.

###
https://www.msys2.org/
###


###

Access Terminal Properties:

Right-click on the MSYS2 terminal title bar.
Select "Options" from the context menu.
Modify Terminal Colors:

In the Options window, navigate to the "Terminal" section.
Click on the "Colors" tab to customize the color settings.
Set Desert Theme Colors:

You need to manually set the colors to match the desert theme. Here are the typical color values for the desert theme:

Background Color: #3F3F3F (R: 63, G: 63, B: 63)

Foreground Color: #DCDCCC (R: 220, G: 220, B: 204)

###


# Update MSYS2 package database
pacman -Syu

# Install the MinGW-w64 GCC Fortran compiler
pacman -S mingw-w64-x86_64-gcc-fortran


nano ~/.bashrc


export PATH="/mingw64/bin:$PATH"
source ~/.bashrc
gfortran --version


