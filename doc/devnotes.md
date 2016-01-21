#tips

##file conventions

- all files should contain one and only one function, unless it is a function that is **only ever** called from another function in that file.
- all files should be have the same name as the main function it contains.



##function naming

all exported functions should have the 'bioarch_' prefix so that it is easy to find help


##required packages

to make sure a required package gets installed when you install a package, all you have to do is use

      require(package)
  
intead of 

      library(package)
