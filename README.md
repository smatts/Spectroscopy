[![CC BY-SA 4.0][cc-by-sa-shield]][cc-by-sa]    [![DOI](https://zenodo.org/badge/390316709.svg)](https://zenodo.org/badge/latestdoi/390316709)

This work is licensed under a
[Creative Commons Attribution-ShareAlike 4.0 International License][cc-by-sa].

[cc-by-sa]: http://creativecommons.org/licenses/by-sa/4.0/
[cc-by-sa-image]: https://licensebuttons.net/l/by-sa/4.0/88x31.png
[cc-by-sa-shield]: https://img.shields.io/badge/License-CC%20BY--SA%204.0-lightgrey.svg


# Experiments in Spectroscopy

[These](https://github.com/MarkusLippitz/Spectroscopy/blob/to-be-published/Spectroscopy.pdf) are the lecture notes of my lecture on optical spectroscopy. The lecture aims at students in the first year of the master's program, but should be accessible also to students in the last year of the bachelor's program. You need some quantum mechanics and an introduction to the physics of molecules.

The idea of this lecture was to focus on _experiments_. Each chapter is build around an experiment and tries to explain everything that is needed to appreciate that experiment. As I am convinced that one only understands something if one really has worked with the concepts, each chapter has a _task_, which is in many cases to evaluate 'real' experimental data, or to simulate the experiment in a computer. The chapters are rather independent of each other, although the order has some sense. One thus can choose appealing chapters and skip others.

In the corona semester, I handed out these lecture notes to the students and prepared short (approx. 15 minutes) introductory videos to each chapter. I have the impression that it is easier for me to put things into context when talking then when writing, as I dare to be more sloppy. I did not enforce any timing, any rhythm during the semester. Everyone could select his or her pace and choice of chapters. During live video conferences, we discussed open questions and the results of the students work on the tasks. The tasks were also uploaded on the elearning platform and visible to all others who have completed a task.

These notes are 'work in progress', and probably never really finished. If you find mistakes, please tell me. I am also always interested in other sources covering these topics.
The most current version of the lecture notes can be found at https://github.com/MarkusLippitz/Spectrosocpy. I have put everything under a CC-BY-SA license. In my words: feel free to do with it whatever you like. If you make your work available to the public, mention me and use a similar license. 


The lecture notes are typeset using the LaTeX class ['tufte-book'](https://tufte-latex.github.io/tufte-latex/) by Bil Kleb, Bill Wood, and Kevin Godby, which approximates the work of [Edward Tufte](https://www.edwardtufte.com/). I applied many of the modifications introduced by Dirk Eddelbuettel in the ['tint'](https://dirk.eddelbuettel.com/code/tint.html}) R package. For the time being, the source is LaTeX, not markdown.

## Structure of the repository

The chapters are individual tex-files in the directory tree and included in the main tex file. Figures are mostly made by tikz. Tikz-external is used to save compilation time. It generates pdfs in the figures tree. Some chapters contain handout material for exercises, which is in the handout subdirectory of the chapter.

