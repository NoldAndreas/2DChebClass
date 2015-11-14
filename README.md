# README #

2DChebClass is a code Andreas Nold and Ben Goddard developed as a project for some simple DFT computations in 2012 by at Imperial College London. Since then, it evolved into a library of classes and functions to solve 1D and 2D DFT and DDFT problems. 

## How to get started

After downloading the version control system 'git', run, in the folder you want your git repository to be created
```
$ git clone https://NoldAndreas@bitbucket.org/NoldAndreas/2dchebclass.git/wiki
```

This should download the content of the code to your computer. 

The first file to edit once you have downloaded the code is "AddPaths.m" in the main folder. 
Please add an option to the switch statement to identify your computer via its MAC address. Also, define via "dirData" a folder where the computational results should be saved.

```
case '67-CF-65-55-C1-82'  %YOUR MAC ADDRESS
            dirData    = 'D:\2DChebData';    %Location of the data folder.
```

The content in the dirDDFT folder will be uploaded, so don't add data files here, as this means that the data limit we have on bitbucket would be exceeded quickly.

### How do I run tests? ###

* Go to 
* Deployment instructions

## Matlab-Versions
The code currently runs using Matlab2014a.

### What is this repository for? ###

* Quick summary
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### Contribution guidelines ###

 (i.e. 'don't merge onto the master branch, put your name on the branches, a link to git commands') 

* only edit/ create branches which start with your first name
* merging of personal branches with the master branch is only done by Andreas and Ben 
* computations are performed 

### Who do I talk to? ###
