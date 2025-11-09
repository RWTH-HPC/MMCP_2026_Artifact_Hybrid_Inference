# Getting Started (Useful knowledge)  # {#ugGettingStartedUsefulKnowledge}

[TOC]

The goal of this section is to provide more background information, which is helpful when you work with @maia.

## Prerequisites

Before trying to work with @maia, you should:

* Have some basic **Unix** knowledge, especially in using a command
  line interface (CLI) such as a terminal or konsole with a bash
  shell, just type `Ctrl + Alt + T` in a KDE desktop to open it. If
  you are using a different desktop environment, please read the
  [Survival Guide for UNIX
  beginners](http://matt.might.net/articles/basic-unix/) and this
  [Introduction to
  Unix](http://www.doc.ic.ac.uk/~wjk/UnixIntro/). Besides the basics,
  you should know how to use unix commands such as
  [**grep**](https://www.man7.org/linux/man-pages/man1/grep.1.html)
  with regular expressions. It is also useful to know [how to
  configure your shell](http://www.hypexr.org/bash_tutorial.php).

* Know how to **edit text files** under unix with your preferred editor. The best
and most common editors are **emacs**, **vim** and **kate** (some PhD
students at AIA may have a different opinion). Choose
one of them or find an alternative one that suits you even better and
learn the basics such as text searching, copy and paste, adding macros, etc.. For more information, have a look at
  * **emacs**: [A guided tour of
emacs](https://www.gnu.org/software/emacs/tour) and an emacs [cheat
sheet](https://www.gnu.org/software/emacs/refcards/pdf/refcard.pdf)
  * **vim**: [Introduction to vim](http://www.openvim.com) or just type
vimtutor in command line, you might need some time to get used to it,
but this [cheat sheet](https://vim.rtorr.com/) might help.
  * **kate**: [Introduction to
kate](https://docs.kde.org/stable5/en/kate/kate/introduction.html) or
the [The Kate
Handbook](https://docs.kde.org/stable5/en/kate/kate/index.html) might
help.
    
* Know how to write a **computer program** to solve a mathematical problem. The best choice of the programming language in the context of @maia would be **C++**, but other computer languages may also be a good starting point. If you don't have any prior knowledge, the MIT course [Introduction to Computer Science and Programming](http://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-00sc-introduction-to-computer-science-and-programming-spring-2011/) may be helpful for you.

* Learn details of the object oriented programming language **C++**, which is extremely powerful, but also sometimes difficult to understand. You should be able to read existing program lines and know how to write basic C++ code lines. There are hundreds of tutorials available on the world wide web and book shops. A short list of tutorials is here:
    * http://www.cplusplus.com/doc/tutorial/
    * http://www.learncpp.com/
    * https://www.cprogramming.com/tutorial.html
    * Book about what **everybody needs to know about C++** titled [Accelerated C++: Practical Programming by Example](http://159.69.3.96/ebooks/IT/PROGRAMMING/Cpp/Addison_Wesley_Accelerated_Cpp.pdf) by Andrew Koenig, which is over 300 pages. It is short, to the point, well explained, and fast to read **if you know how to program*. 
    * To play around online you can go to http://cpp.sh
    
* Have some basic knowledge of the programming language **Python**:
    * https://www.python.org/about/gettingstarted/
    * https://www.w3schools.com/python/python_getstarted.asp
    
* If you plan to develop parallel algorithms, please have a look at **parallel programming using MPI** (Message Passing Interface), and understand at least the most basics MPI operations, which are explained, e.g., in this tutorial:
    * http://mpitutorial.com/tutorials/
    
* Have some basic knowledge using
     * [Git/GitLab](https://docs.gitlab.com/ee/tutorials/make_your_first_git_commit.html), which is a distributed version control system for collaborative software development. 
     * [SSH](https://linuxhandbook.com/ssh-basics/) (Secure Shell) used to connect to remote computers and transferring data.
     * the Slurm Workload Manager [Quick Start User Guide](https://slurm.schedmd.com/quickstart.html), with some additional commands such as `shosts` and `si` only available at AIA.
     * [HDFView](https://docs.hdfgroup.org/archive/support/products/java/hdfview/UsersGuide/ug02start.html) for opening hdf5-files, which are used for storing the simulation results.
     * [ncdump](https://docs.unidata.ucar.edu/netcdf-c/current/netcdf_documentation.html) for opening netcdf-files, which are used for storing the simulation results.
     * [ParaView](https://www.paraview.org/tutorials/) for visualization and postprocessing, which can be opened with the command `parav`or `paraview`.


A lot of your work will be learning by doing, but some basic knowledge on the aforementioned topics will speed up you work significantly and you won't be totally lost once you start working with or on the code.  

@note
If you are using a CLI in a Terminal with the [**bash**](https://www.gnu.org/software/bash/manual/bash.pdf), you can complete a particular command or filename using the `Tab' key if the name is already uniquely specified. Pressing the `Tab' key twice will give you the list of options available.  You can also use the arrow keys to reuse and edit commands in your command history.

Besides programming and computer knowledge you should also have a good understanding of the physics related to your simulation problem such as **general fluid mechanics**, **acoustics**, **heat transfer** etc.. While this is certainly not a knowledge you can pick up on the fly, there are a lot of text books, papers and internet resources that may help to refresh and improve your knowledge.

* NASA's [Compressible Aerodynamics Index](http://www.grc.nasa.gov/WWW/k-12/airplane/shortc.html) is a collection of articles on various fluid mechanics-related topics.
* [An introduction to Aeroacoustics](https://www.researchgate.net/publication/228690041_An_introduction_to_aeroacoustics)
* [Heat transfer introduction](https://www.researchgate.net/publication/323144799_Heat_transfer_introduction)

