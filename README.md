# ANSYS-USERMAT-CMSG
The conventional theory of mechanism-based strain gradient plasticity is realized via ANSYS user programmable features for USERMAT subroutine.   
Strain gradient plasticity is a theory of continuum solid mechanics which aims at modeling the irreversible mechanical behavior of materials, with specific focus on metals and on their response at appropriately small size, typically on the order of micrometers or less. For small scale, a variation in size leads to a peculiar effect, denoted as “smaller being stronger.” At the crack tip, taking into account gradient plasticity leads to a significant increase of stresses.

<center>

<img src="Images\Crack tip stress.png" width="600" alt="Crack tip stress"> 

</center>

<br>

 In source files you can find a APDL script example for 2D cracked body. After compiling and attaching present dynamic link library as ANSYS user material copy and run [this file](https://github.com/Andrey-Fog/ANSYS-USERMAT-CMSG/blob/main/Source/APDL-%202D%20crack%20example.txt) from ANSYS Mechanical command line.

<br>

 ## Research results  
- [Crack tip fields and fracture resistance parameters based on strain gradient plasticity](https://doi.org/10.1016/j.ijsolstr.2020.10.015)  

- [Mode I and mode II stress intensity factors and dislocation density behaviour in strain gradient plasticity](https://doi.org/10.1016/j.tafmec.2021.103128)

- [Inversion of dislocation densities under mixed mode fracture](https://doi.org/10.1016/j.engfailanal.2022.106311)

## Acknowledgment
I want to say thank a lot to ​​​[​Emilio Martínez Pañeda](https://www.empaneda.com/) for the detailed presentation of his methods in open acsess and active assistance. This helped a lot in the implementation of CMSG theory into ANSYS. This project is reworked his ABAQUS subroutine presented [here](https://doi.org/10.1016/j.ijsolstr.2015.02.010). 


<br>
<br>
<br>

## Instructions for compiling and attaching USERMATLIB.DLL 

---

Install Visual Studio first and then Intel fortran compiler. When installing the compiler, select "Integrate into Visual Studio". Supported versions can be found in the ANSYS documentation in the section on User Programmable Features (UPF). Add LIB and INCLUDE variables in the system environment. Create new solution and add new fortran dll project. The name of the created library must be "USERMATLIB.DLL". Add all fortran files from Source directory to your dll project. Tune compiler according to instructions present below. After compiling connect library to ANSYS.


<br>

### Connecting to ANSYS

---

After creation the dll file you have to connect this library to ANSYS:

<br>

**1. Create environment variable named ANS_USER_PATH**

*My Computer->Properties->Advanced system settings->Advanced*  

On the tab, click on the button:

*Environment Variables->System Variables->New*

<br>

**2. In the variable value field, specify the path to the folder where library is located. Use only latin characters in the path.**

*For example:* 
>C:\Username\......\Usermatlib
   
If everything is connected correctly in the ANSYS output window at startup there will be a line 

```
User link path <ANS_USER_PATH>: *path to your folder*" 
```
<br>

**3. After launching the ANSYS, create an user material**

*Preprocessor->Material Props->Material models*

<br>

**4. In the drop-down list of materials, select**

*Structural->Specialized Materials->User material options->user material*


And add cells. There should be 6 properties in total. Of which:

| NN  |     | Property                              |
| --- | --- | ------------------------------------  |
|  C1 |  -  |Young modulus                          |
|  C2 |  -  |Puasson ratio                          |
|  C3 |  -  |Yelding stress                         | 
|  C4 |  -  |Intrinsic material length              |  
|  C5 |  -  |Strain hardening exponent (0 < N < 1)  |
|  C6 |  -  |flag, 0 or 1. [Usually](https://doi.org/10.1016/S1359-6454(99)00020-8) 1 for metals|

In command line it will be looks like present bellow

>!* Define parameters related to CMSG model  
>!* Modulus of Elasticity  
>Young	= 200000   
>!* Poisson ratio  
>nu	= 0.3  
>!* Yield Strength  
>S02	= 300  
>!* Length parameter in the CMSG model (approximate order 1e-6 meters)  
>Leng 	= 1E-6  
>!* Strain hardering exponent (0 < N < 1)  
>SHE	= 0.15  
>!* Flag  0 or 1. Usually 1 for metals  
>Flag 	= 1
>     
>!* add user model  
>TB,USER,1,1,6,  
>TBTEMP,0  
>TBDATA,,Young,nu,S02,Leng,SHE,Flag  

**5. Add 13 state variables**  

*Preprocessor->Material Props->Material models->Structural->Specialized Materials->User material options->State Variables*

| SVAR| Value                                 |
| --- | ------------------------------------- |
| 1-3 | Coordinates of integration points     |
| 4-9 | Derivatives x,y,z,xy,yz,zx            |
|  10 | Gradient value                        |
|  11 | Statistically necessary dislocations  |
|  12 | Geometrically necessary dislocations  |
|  13 | Total density of dislocations         |

APDL script for preprocessor section

>TB,STAT,1,1,13,  
>TBTEMP,0  
>TBDATA,,0,0,0,0,0,0  
>TBDATA,,0,0,0,0,0,0  
>TBDATA,,0  

**6. Access to user arrays**

Before starting on the solution in the solver (/SOL) in the command line of ANSYS, write the line:

- to save every substeps results  
> OUTRES,SVAR,ALL

- to save only the last step  
> OUTRES,SVAR,LAST

In order for all elements of user arrays to be available, command GRAF must be used in the postprocessor (/POST) section.  

> /GRA,FULL

That's all. Further we work as with the usual scheme.

<br>
<br>
<br>

#### COMPILATOR SETTINGS (*projectname*->properties). 

| Name     |   | Value |
| ----------- | ----------- |----------- |
|Supress startup banner:| 	 - |	            Yes (/nologo) | 
|Additional include Directories:|  - |	        C:\Program Files\ANSYS Inc\v***\ansys\customize\include | 
|Optimization:| 			 - |	            Disable (/Od) |
|Preprocessor definitions:| 	 - |	        /DNOSTDCALL /DARGTRAIL /DPCWIN64_SYS /DPCWINX64_SYS /DPCWINNT_SYS /DCADOE_ANSYS /D__EFL /DFORTRAN /auto /c /Fo.\ /MD /W0  |
|Debug information Format:|		 - |            Full (/debug:full)  |
|Preprocess Source file:|		 - |            Yes (/fpp)  |
|Preprocessor Definitions to fpp only:| -|	Yes (/noD)  |
|Use Portlib Library:| 		 - |	            Yes (/4Yportlib)  |

#### LINKER SETTINGS  
| Name    |  |    Value |
| ----------- | ----------- |----------- |
|Enable incremental linking:| - |		No (/INCREMENTAL:NO)  |
|Supress startup banner:|  - |		    Yes (/nologo) | 
|Additional library Directories:|  - |	C:\Program Files\ANSYS Inc\v***\ansys\custom\lib\winx64|  
|Additional dependencies:| 	 - |	    ANSYS.LIB  |
|Generate debug info: |	 - |		    Yes (/DEBUG)  |

*** - your version of ANSYS
All other settings by default. Its allows me connect to ANSYS for debugging.
