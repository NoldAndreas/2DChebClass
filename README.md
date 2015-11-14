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

## Matlab-Versions
The code currently runs using Matlab2014a.

### What is this repository for? ###

* Quick summary
* Version
* [Learn Markdown](https://bitbucket.org/tutorials/markdowndemo)

### How do I get set up? ###

* Summary of set up
* Configuration
* Dependencies
* Database configuration
* How to run tests
* Deployment instructions

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact