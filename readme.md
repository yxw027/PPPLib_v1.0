# Precise Point Positioning Library
![N1hvPe.png](https://s1.ax1x.com/2020/06/21/N1hvPe.png)
## 0. 更新日志
    2020/06/21 PPPLib v1.0发布，支持GBERJ单频-三频非组合/组合PPP解算

## 1. 关于PPPLib
Precise Point Positioning Libary(PPPLib)是一款采用C/C++语言编写的多模多频GNSS-PPP数据处理算法库；主要功能包括：
	
	× GNSS常用文件读取
	+ 支持GPS,BDS,GAL,GLO,QZS
	* 任意单系统或任意系统间组合
	× 标准单点定位SPP
	+ 精密单点定位PPP
	* 多系统单频/双频/三频 非组合PPP
	+ 多系统单频/双频/三频 无电离层组合PPP
	- 电离层约束非组合PPP
	× 支持BDS-3卫星及其新信号
	* 丰富的结果输出
	+ 丰富的调试信息输出 
    × 数据下载/结果分析
PPPLib定义为入门级GNSS-PPP数据处理软件，是作者本人在学习GNSS-PPP数据处理过程中的程序化实现，因此程序的计算性能和误差处理模型仍有较大提升空间。
但整体上实现了多模多频PPP数据处理，希望对刚入门的PPPer有所帮助。目前PPPLib正在进行PPP/INS融合算法的编写，后续会不定期更新，规划的开发路线为:	
    
    0. PPP/INS多源融合(松紧)算法
    1. PPP-AR算法
    2. PPPLib-GUI开发
    
## 2. 如何使用
PPPLib采用Cmake进行工程管理，目前软件仅支持在linux下的编译和运行，Windows暂未做测试。测试系统包括虚拟机Ubantu16.0以及Ubantu16.0系统，
作者在开发过程中使用的是JetBrains CLion2019.3软件平台，尚未做其他软件平台的测试。(依作者经验，在Linux下做开发比较好，遇到问题更容易解决，所以推荐使用Linux)。
1. 克隆：打开终端 
    cd [Your Path] \
    git clone https://github.com/heiwa0519/PPPLib.git
2. 编译： \
    cd [Your Path]/PPPLib \
    mkdir build \
    cd build \
    cmake .. \
    make 
3. 运行： \
    cd [Your Path]/PPPLib/bin \
    ./PPPMain -pd 2019/12/01 -do 1 -level 128 -sys GBERJ -md PPP-KINE -ion 4 -frq 1   \
    具体配置信息可以查看用户文档([Your Path]/PPPLib/doc/PPPLib用户手册) 
一组示例数据上传至百度云，扫码下载后放至项目根目录下即可

[<div align=center>![U9KKHO.th.jpg](https://s1.ax1x.com/2020/07/05/U9KKHO.th.jpg)](https://imgchr.com/i/U9KKHO)

## 3. 待完善
总结了目前PPPLib有待改进的点，个人能力有限，如果有好的建议可以直接联系本人，一起完善PPPLib。 
    
    0. C++11代码编写不规范
    1. 更加稳健鲁棒的PPP算法
    2. 多系统非组合计算效率低
    3. 完善电离层约束PPP
    4. 完善GLONASS伪距频间偏差IFCB估计
    5. 鉴于BDS全球组网在即，BDS的误差模型仍需精化和完善
    6. 实测动态数据测试(如有能提供实测GNSS数据(长时间),GNSS/INS/里程计数据的请联系作者,非常感谢)   
        
## 4. 致谢
首先致敬[RTKLIB](https://github.com/tomojitakasu/RTKLIB/tree/rtklib_2.4.3)定位软件作者Tomoji Tkasu（高须知二）先生，其无私的开源精神和优雅的程序设计令我敬佩不已。
PPPLib的实现同时还参考了山东大学Mowen Li的[HPRTK](https://github.com/Bemo12)软件,
山东科技大学周锋老师的[GAMP](https://link.springer.com/article/10.1007/s10291-018-0699-9)软件，测地所肖恭伟博士的[MGAPP](https://github.com/XiaoGongWei/MG_APP)软件以及西北工业大学严恭敏老师的PINS软件，PPPLib的日志输出使用的是
[easyloggingpp](https://github.com/amrayn/easyloggingpp)再此一并感谢。同时感谢师弟朱霆在文档撰写方面的协助。开发不易，如有引用，敬请注明出处。 

## 联系作者
在使用过程中有任何问题或者好的建议请联系作者；PPPLib将一直维护下去，但个人精力和能力有限，PPPLib-dev-Group欢迎感兴趣的同学加入。\
地址： 中国矿业大学环境与测绘学院 陈超 \
QQ: 565681993 \
e-mail: cchen@cumt.edu.cn

**Supported by [Guobin Chang's Lab](https://www.researchgate.net/lab/Guobin-Chang-Lab)**

****************************************************************
## 0. Revision Logs
    2020/06/21 PPPLib v1.0 released, supported GBERJ single- to triple- frequency ionosphere-free/uncombined PPP
    

## 1. About PPPLib
The PPPLib adopts post-processing mode for PPP processing to integrate multi-frequency and multi-
GNSS data (GPS, BDS, Galileo and GLONASS as well as QZSS). PPPLib was developed in C/C++
environment, it is recommended that compiles and runs it in Linux system. The main features include:
    
    * GNSS common functions
    + Support GPS/BDS/GAL/GLO/QZS
    - Support any single system or combined
    × Standard single point positioning (SPP)
    * Precise point positioning (PPP)
    + single- to triple frequency with multi-constellation PPP using ionosphere-free observations
    - single- to triple frequency with multi-constellation PPP using uncombined observations
    x inospheric constranit PPP
    + Support BDS-3 satellites and new signals
    * Ample output
    - Ample debug information
    + data download and result analysis

Future features of PPPLib will be:
    
    0. loosely/tightly coupled PPP/INS (on the way)
    1. PPP-AR
    2. PPPLib-GUI
 
## 2. How to use
PPPLib uses Cmake for project management. Currently the software supports compile and run in Linux and Windows. 
The test system includes the virtual machine Ubantu 16.0, Ubantu 16.0 and Win10. The author uses JetBrains CLion2019.3 software platform for developing and debug. No test other platforms.
    
    open the terminal 
    0. Clone: 
       cd [Your Local Path]
       git clone https://github.com/heiwa0519/PPPLib.git 
    1. Compile:
       cd [Your Local Path]/PPPLib
       mkdir build
       cd build
       cmake..
       make
    2. Run:
       cd [Your Local Path]/PPPLib/bin
       ./PPPMain -pd 2019/12/01 -do 1 -level 128 -sys GBERJ -md PPP-KINE -ion 4 -frq 3
For more specific configuration information, please refer to PPPLib user manual.     

Upload a set of sample data to Baidu Cloud, please scan the code to download. Then place the example folder in the root directory

[<div align=center>![U9KKHO.th.jpg](https://s1.ax1x.com/2020/07/05/U9KKHO.th.jpg)](https://imgchr.com/i/U9KKHO)


## 3. To be improved
    
    0. No standard C++ coding
    1. Need more robust PPP algorithm
    2. Low wfficiency for multi-GNSS and multi-frequency uncombined PPP
    3. Bug in ionospheric constraint PPP
    4. Bug in GLONASS IFB estimation
    5. In view of the imminent BDS global networking, the error modeling for BDS satellites should be imporved
    6. Lack dynamic test (if your can provide any measured GNSS data (long time), or GNSS/INS/odometer data, please contact the author, thank you very much)

## 4. Acknowledgement
First of all, I pay tribute to Mr. Tomoji Tkasu, the author of [RTKLIB](https://github.com/tomojitakasu/RTKLIB/tree/rtklib_2.4.3) software. I admire him for his selfless
open source spirit and elegant programming. The developing of PPPLib also refers to [HPRTK](https://github.com/Bemo12) of Mowen Li from Shandong University,
the [GAMP](https://link.springer.com/article/10.1007/s10291-018-0699-9) software of Feng Zhou from Shandong University of Science and Technology,
the [MGAPP](https://github.com/XiaoGongWei/MG_APP) of Gongwei Xiao from Institute of Geodesy and Geophysics，Chinese Academy of Sciences, and the
PINS software of Gongmin Yan from Northwestern Polytechnic University. The log system of PPPLib uses [easyloggingpp](https://github.com/amrayn/easyloggingpp).
Thanks to the authors of the above software. It is not easy to develop, please indicate the source if you have any references.

## Contact Author
Any suggestions, corrections, and comments about PPPLib are sincerely welcomed and could be sent to: \
Author: Chao Chen \
QQ: 565681993 \
E-mail: cchen@cumt.edu.cn \
Addres: School of Environment and Geo-informatics, China University of Mining and Technology

**Supported by [Guobin Chang's Lab](https://www.researchgate.net/lab/Guobin-Chang-Lab)**   
       
    


