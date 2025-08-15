
# Rosalind Problems

[Rosalind.info] is a site for learning bioinformatics at one's own pace. It contains a bunch of programming questions to practice writing algorithms for bioinformatics. This repo contains my solutions to these problems. I got hooked on these problems in summer 2025 after doing the first 16 within 24 hours of finding the site.



### Files

Since the problems are arranged into a [hierarchy of difficulty](https://rosalind.info/problems/tree-view/), I've split the solutions into files according to what "level" of the problem tree they appear in. Note that the first three levels contain just one problem each, since these are the introductory problems.

- `LevelN.R`
	- These are the main solutions files. Each problem has a title, a function to solve it, and a call to that function for convenience.
	- These files also contain one-off helper functions. These appear just before the problems they are utilized in.
- `helper_functions.R`
	- This file contains helper functions used in multiple problems, for example a `fasta` parser.
- `LevelN_XXXX.py`
	- These files contain Python solutions to select problems (see FAQ). Each file is the solution to one problem, unlike the R files above. The filename includes both the "level" of the problem and the 2-4 letter problem code used by Rosalind.
- `useful_charts`
	- This folder contains some charts provided by Rosalind to help solve the problems.

### FAQ

- **Why did you use R *and* Python?**
	- I used R for the majority of the problems. I'm a statistician by training, and R is the language I'm most comfortable with by far. That being said, R is not without its downsides. If a problem required an object-oriented approach for example, I used Python instead.
- **Rosalind has a lot of problems on it. Which ones are these?**
	- These problems are from the "Bioinformatics Stronghold" problem bank. As previously noted, the problems are arranged hierarchically, so you cannot attempt a problem until you've successfully completed it's prerequisites. A visualization for this hierarchy is [here](https://rosalind.info/problems/tree-view/), or you can see the problems in list view [here](https://rosalind.info/problems/list-view/).
- **Aren't there solutions to these problems on the internet already?**
	- Yes. You can find plenty of others' solutions out there if you wish. I'm putting mine out there purely as a demonstration of my own skills.
- **Are these solutions entirely your own?**
	- Yes. I did not consult other people's solutions when writing my own. I permitted myself reasonable use of outside resources such as StackExchange and AI tools (for example, to help optimize my solutions to fit within the 5-minute time limit) but I did not use any outside resource that would trivialize the problem-solving process. Sometimes I looked at other folks' solutions *afterward*, to see how other people approached it. So, I am very much aware that there are better ways to solve these problems, but I've tried to preserve the original state of my solutions as much as possible. 