
#Basics of R for non-programmers

##SUMMARY

1. Running R / RStudio
- Working with packages
- learning R with swirl
- Installing github packages
- the bioarch package
- writing packages with devtools
- tips

---

##RUNNING R

_Assuming R is installed..._

First thing, we need to have a look at the different ways of running R

- from a command line
- in Rstudio

RStudio is a great way to start, especially if you are used to GUIs. But there are occasions (such as when you are working remotely) where the command line is useful. 

OK, let's start RStudio. For more help, see:

[https://www.rstudio.com/resources/training/online-learning/]

---
##WORKING WITH PACKAGES




There are five levels of programming in R:

1. **your own code**		Programs you'd only use yourself, for analysing and presenting your own data. Before we do _any_ programming, we should generally check if a solution exists elsewhere!
- **shared code (bad practice)**	Programs you share locally with colleagues - proprietary stuff. It is very easy to email programs to one another, but often this can cause problems if the programs don't work for one reson or another. It is very easy to lose track of the 'correct' version of your code, and it is often poorly documented. 
- **packages**	It's not much more work to group your programs together into a package that you develop as a group. We'll show how to do this using github below. The advantages are that it is easier to document what you are doing as you go along, and you can timestamp versions of the code so it is easy to keep track of how things work at any point in the development of the software. 
- **CRAN packages**		Next level ('The Comprehensive R Archive Network'). Specialist programs, maintained by the R community. Rigorously tested, very powerful, but sometimes not very well documented and quirky to use. The requirements for CRAN hosting can be quite challenging to meet, but that's why it's the 'gold standard'!
- **R 'Core'**		This is the R programming language, operators functions etc. We'll have a look at how to learn this stuff over the next few sections


###INSTALLING PACKAGES

Most of the time, you'll be installing packages from CRAN, since these are the most tried and tested ones. These are the easiest packages to install. 

I'm going to show you two packages that you should know about:

- **Swirl** - Friendly R tutorials
- **devtools** - ways to help you oragnise your code into packages from the outset. (More of this)

---
##Learning R with Swirl

Packages can be about anything. 'Swirl' is a good package to start with, because it contains tutorials for learning R.

Let's install swirl using the RStudio menu: 
   *tools->install packages...*
gets you to a dialog box. Type ```swirl``` into the 'Packages' text-entry field. 

Alternatively, at the command line (in the bottom left panel of RStudio) type

```
   install.packages("swirl")
```

**REMEMBER!** *Installing* a package is not the same as *loading* a package so that you can use it. Think of installing a package as like buying a book and putting it on your bookshelf. You can't read it while it's on the shelf - you have to hold it first. Similarly, you have to load the package before you can use it, using the `library` command:

```
library("swirl")
```

For more on swirl, see:

[http://swirlstats.com/students.html]

This is a set of interactive tutorials. Let's install that and go through a bit

There's also an online course starting 18th Jan 2015 - it's free - coursera are pretty good:

[https://www.coursera.org/learn/r-programming]



---
##Installing github packages

Packages that are not hosted by CRAN require a little more work to install. Fortunately, the 'devtools' package helps with this. Devtools *is* hosted by CRAN, so we can install it in the same way we installed Swirl. 

devtools is a great thing for building packages. If you plan to share your software, you should **begin** by creating a package. It will save you time in the long run! If you've ever come back to some code, or even a complicated spreadsheet and wondered how it works, you should know what I mean. Here's a great blog post about writing packages (if a little out of date):

[http://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/]

To demonstrate this, we're going to load the `bioarch` package that this tutorial came with.
Here's the command-line for installing and loading devtools (you can use the RStudio menu if you like):
 
 `> install.packages("devtools")`
 
 `> library(devtools)`
 
Then you can use the `install_github` function to install bioarch:
 
 `> install_github("franticspider/bioarch")`
 
Of course, you'll only need to do that once for your R installation. After that, whenever you want to use bioarch, simply enter: 
 
 `> library("bioarch")`

We'll talk further about this next time, but rember, if a package is written entirely in R, then install will work no matter what platform you are using. Extra steps have to be taken if the R package uses programs written in other languages (yes, R can do that)
For example, the q2e package we've built uses Julie Wilson's C code, so we have to take a few more steps to make it work in Windows


---
##Tips

###LISTING PACKAGES

installed.packages() 

will list everything.

###REMOVING PACKAGES.

should be simple:

remove.packages("devtools")

If the package isn't there, you'll get this error:

	Error in remove.packages : there is no package called ‘devtools’

sometimes the package to be removed isn't in the usual directory. If you think it's still there do this:

.libPaths()

which will show all the paths that packages *can* be on. You can list these directories in the remove step like this:

remove.packages("devtools","/usr/lib/R/library")




