
We will be working with bash. Bash is a command line language, meaning that you type a command to do something (usually followed by some arguments), hit return/enter. Bash reads this command (and what you have added next to it) and runs it. Some examples that are simple commands are:
- open and print file X
- print a list of the files that are in directory X
- print the working directory
- sort the contents of file X
- copy/move file X to position Y
Refer to the powerpoint presentation for some examples and how to navigate yourself in the linux environment.

The exercises below aim to familiarize yourself with linux commands. See linux commands as the tool to extract the biological information from the file. We will start with simple commands and gradually move to more complicated ones. The goal is not to get the answer but to think of which commands you should be using, with what sequence and with what arguments. Feel free to go to the manual pages of each command, try different syntaxes and different arguments. You will not ruin your computer and don't get afraid if a very big output floods your screen.

If you get stuck, you can see the answer of the question but it is important to read it and understand what's happening and why it is written that way. Again, the goal of these exercises is not to get to the answer but to understand how we get there. For pipes, to understand how data are processed sequentially, it is always helpful print the output of each step, see what's happening and then write the next command(s).

We will be working with the file, downloaded from BioMart: mart_export.txt.gz

Notes:

_i)_ For what we will cover here, you will need to consult the 'man' (manual) pages. You can do this by typing
  ```
  man wc
  man sort
```

etc. and then hitting enter/return. The manual pages contain useful examples. Once in the manual page, type 'q' to quit and return to the command line.

_ii)_ The general syntax in bash is:
```
<command> <argument(s)> <file(s)>
```

The command can be something like 'wc' (for Word Count), argument(s) are used to tell the command to do something specific (e.g. count only the lines) or to specify something (e.g. the desired column). So, let's get right to it with our first tasks:


# **EXERCISE 1**

Download the file:
wget https://github.com/ceccarellilab/CancerBioinformaticsCourse/raw/refs/heads/main/Biological%20DB%20\(Aris\)/mart_export.txt.gz

1) Uncompress the file (hint: gunzip)

2) Get the number of lines, words and characters in the uncompressed file (hint: wc)

3) Get only the number of lines in the file (hint: wc -l )

4) Get the first 10 lines of the file (hint: head)

5) Get the last 10 lines of the file (hint: tail)

6) Extract the second column (field) from the file (hint: cut)

7) Extract the lines that contain CLOCK (hint: grep)

8) Extract the lines that contain the word TP53 (hint: you have to use a specific argument with grep)

_Great! So far, you've made significant progress in bash! Let's make some more. With the commands above, you were seeing the output on your screen. Linux offers a very simple way to work on that output. For example, we can get all the instances of CLOCK, extract the transcript types (6th column), extract only the protein coding instances and count them. You know the commands for all these steps and now we are going to learn how to put these commands in order: linux uses pipes, the symbol '|'. This symbol tells linux to take the output of the previous command and use it as the input of the next command (so you don't have to specify a file). The syntax is:_

```
 grep CLOCK mart_export.txt | cut -f 7
```

_Try that and see what you get._

9) Get the 135th line of the file
_Hints: you need to use a combination of 'head' and 'tail'. You can use multiple pipes, one after the other, e.g.:_
	<command> file | <other_command> | <yet_another_command> | <one_more_command>
_How to do this: Get the first 135 lines and then get the last one_

Now, you can complete the task:

10) Get all the instances of CLOCK, extract the transcript types (5th column), extract only the protein coding instances and count them. 
Tip: instead of "grep CLOCK mart_export.txt" you can use "cat mart_export.txt | grep CLOCK"

11) Get all the instances of CLOCK, extract the transcript IDs and sort them (hint: sort)

12) Get all the instances of CLOCK, extract the transcript types, sort them and filter out repeated lines (hint: uniq)

13) Get all the instances of CLOCK, extract the transcript types and filter out repeated lines.

_If you followed the instructions of 12 and 13 exactly word-by-word, then the outputs should differ. Specifically, 13 should have returned three instances of "protein_coding", three instances of "retained_intron" and one instance of "processed_transcript". Why is that the case? (hint: look at the manual of uniq)_

_Now, by getting the appropriate columns, sorting, filtering out repeated entries and counting, you should be able to answer the following questions:_

14) How many different transcript types does TP53 gene encodes for? Can you count how many different transcripts are of each type? (hint: uniq -c). Can you sort the output numerically? (hint: sort -n)

**_Time for some harder questions:_**


15) How many protein_coding genes are in the human genome? (hint: don't confuse gene type with transcript type and don't double-count genes with multiple transcripts, the answer is 20096)

16) How many miRNA genes are in the human genome? 

17) Which chromosome has the highest number of protein coding genes? (Answer: chromosome 1 with 2066 coding genes)

18) How many distinct gene names (column 7) are in the file?
    18.1) How many of them are duplicates? (hint: look at uniq's manual, answer is 21718)
    18.1) How many of them have the letter A?

19) Get the gene IDs (column 1) of all the mitochondrial tRNAs (Mt_tRNA) and save them to a file.

_Hints: build your pipe as you have learned already and to save the output in a file write "> mt_tRNA_IDs.txt" at the end. ">" tells bash to redirect the output not on the screen but on file named "proteinCodingIDs.txt". See the first ten lines and count the number of lines to make sure you got it right._

_Now that we have a file with the IDs of the protein coding genes, we can just use it with a single grep to get them from the file instead of making the pipe each time. That's how it's done (read grep's manual for details):_

```
  cat mart_export.txt | grep -F -f mt_tRNA_IDs.txt
```

20) Open the file with cat and filter the lines to keep the entries that are on chromosome 15

_Hint: awk
I hope you will end up loving awk! It is an elegant and very easy way to filter, do arithmetics and manipulations. Awk reads the file line by line, checks if the specified conditions apply and does something. We specify (i) the conditions (e.g. if chromosome is 15) and (ii) what that 'something' is - it could be to just print the whole line or we could write a whole program. This is what "pattern-action" means in awk's manual. Awk works with fields. So far, we have been using columns as the fields separated by tabs. Awk considers white spaces also as separators. Note that for the two lines below, awk sees two fields, the first one is "CLOCK" and the second is "protein_coding":_
     CLOCK protein_coding
CLOCK          protein_coding

The command below prints only the second and the fifth fields:
```
cat mart_export.txt | grep CLOCK | awk '{print $2,$5}'
```

The command below filters for the protein coding transcripts of CLOCK (the 5th field equals to string "protein_coding")
```
cat mart_export.txt | grep CLOCK | awk '$5=="protein_coding"'
```

And the command below combines the above two awks:
```
cat mart_export.txt | grep CLOCK | awk '$5=="protein_coding" {print $2,$5}'
```

_In programing, "a=2" means to assign the value 2 to a variable stored in memory with the name 'a'. But "a==2" is an expression that checks if the variable is equal to 2 - instead of integers, you can have strings, as above_


_**Time for some more complicated questions**_

21) Filter for the protein coding genes (column 8) on chromosome 16 and count how many of them have transcript types as "retained_intron". (hint: use two awk commands with pipes, answer is 560). Count how many protein coding genes are on chromosome 16.

22) Repeat the above steps for chromosome 20. What do you observe?

_Extra points: make a file with the number of genes per chromosome and another with the number of genes having retained introns as transcript types. Plot these numbers in an excel - what do you observe?_

23) Which protein-coding gene has the most coding splice variants? (Answer/hints: MAPK10 with 116 distinct coding transcript IDs)

24) Which protein-coding gene has the most non-coding splice variants? (Answer: FANCL with 123 distinct transcript IDs that do not code for a protein)

25) How many protein coding genes have coding AND non-coding transcripts? (Hint: if a gene is protein-coding in column 3, then there is at least one protein-coding transcript in column 5 - this means that if there are more than protein-coding transcripts, then the gene ID is repeated in multiple lines)

_**And now back to some simpler questions**_

26) Get the transcript IDs of CLOCK and replace all "0" with "1" (hint: tr)

27) Get the transcript IDs of CLOCK and replace "ENST" with "id" (hint: sed 's/ENST/id/' )

28) Get the transcripts of CLOCK and get their length (hint: do a subtraction and print through awk)
Can you print it like?
 The length of X is Y
where X is the transcript ID and Y is the length

29) How many TP52 transcripts are longer than 10,000 bases? (Answer: 23)

30) Go back to Question 18 and find how many distinct gene names (column 7) are in the file. Then:
    30.1) How many of them have a dot (.)? (hint: google for "escape character" and "how to grep a dot in bash")
    30.2) How many of them contain the letter A followed by any character and then followed by 2? (hint: what does a dot represent? The answer is 1043)


Congratulations! If you made it through all commands, then you are an advanced linux user!



