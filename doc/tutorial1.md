
#Basics of R for non-programmers

##SUMMARY

1. Running R / RStudio
- Installing CRAN packages
- learning R with swirl
- Installing github packages
- the bioarch package
- writing packages with devtools

---

##RUNNING R

_Assuming R is installed..._

First thing, we need to have a look at the different ways of running R

- from a command line
- in Rstudio

RStudio is a great way to start, especially if you are used to GUIs. But there are occasions (such as when you are working remotely) where the command line is useful. 

OK, let's start RStudio

[https://www.rstudio.com/resources/training/online-learning/]

------
##ADDING AND REMOVING PACKAGES

There are five levels of programming in R:

- R Core			This is the R programming language, operators functions etc. 
- CRAN packages		Next level ('The Comprehensive R Archive Network')
- Other packages	Often it's easier to make more specialist packages public outside of CRAN
- local packages	Packages you share on a local network - proprietary stuff
- your own code		Programs you'd only use yourself. 

Before we do _any_ programming, we should generally check if a solution exists elsewhere!



-----
##INSTALLING PACKAGES

Most of the time, you'll be installing packages from CRAN, since these are the most tried and tested ones. 

I'm going to show you two packages that you should know about:

Swirl - Friendly R tutorials
devtools - ways to help you oragnise your code into packages from the outset. 

(We'll come back to swirl later)


------
A really good way to get to know R basics is to use swirl:

[http://swirlstats.com/students.html]

This is a set of interactive tutorials. Let's install that and go through a bit






------
FILES

It gets tiring writing the same stuff over and again. 

Easier to write the code in a file - and load that using

source("file.R")



------
FUNCTIONS

often you want to manipulate data in the same way. this is done with *functions*:

```
addsquare <- function(x){
	y <- x*x
	z <- x+y
	return (z)
}
```

------
DEVTOOLS

devtools is a great thing for building packages. I'm talking about this now because 

http://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/

NB: If a package is written entirely in R, then install will work no matter what platform you are using. 

Extra steps have to be taken if the R package uses programs written in other languages (yes, R can do that)

For example, the q2e package we've built uses Julie Wilson's C code, so we have to 

------
LISTING PACKAGES

installed.packages() 

will list everything.

------
REMOVING PACKAGES.

should be simple:

remove.packages("devtools")

If the package isn't there, you'll get this error:

	Error in remove.packages : there is no package called ‘devtools’

sometimes the package to be removed isn't in the usual directory. If you think it's still there do this:

.libPaths()

which will show all the paths that packages *can* be on. You can list these directories in the remove step like this:

remove.packages("devtools","/usr/lib/R/library")




