# simd-dna
A Python simulator for the SIMD||DNA model of computation, as defined in https://doi.org/10.1007/978-3-030-26807-7_12 (preprint: http://users.ece.utexas.edu/~soloveichik/papers/simd-dna.pdf)

# Setup
1. Install PyCharm: https://www.jetbrains.com/pycharm/
2. In the Welcome screen, select Open (if using PyCharm for the first time) or go to File > Open (if using PyCharm with a pre-existing project loaded) and open the root directory of the repository.
3. Go to File > Settings > Project > Python Interpreter. Choose a Python interpreter from the list if you already have any, or click on the gear icon to the right and choose Add. I recommend using a Virtualenv environment to keep it self-contained.
4. At the bottom, click on Terminal. Type and enter `pip install -r requirements.txt` to install all required packages.

# Usage
## Running
In the project directory, open the context menu for `main.py` and select Run.

## Loading files
To load an existing SIMD||DNA configuration, supply the JSON file name as an argument (either absolute or relative paths are fine.) In PyCharm, go to Run > Edit Configurations... and enter the arguments in the Parameters field. When running for the first time, you may need to set up the configuration first. You can load some of the pre-existing samples such as `rule11.json` or `increment.json` and select Run simulation to see how it works.

## Options
Note: As of this writing, there is no input validation, so any parsing errors or references to non-existent values will throw a runtime exception. This will be addressed in a future update.

The simulator is currently terminal-based, and all options can be accessed by typing the option number in the terminal and hitting Enter.

1. Add cell type  
Input values: Cell name (string), Domain names (comma-separated strings)  
Example:  
`Bit`  
`1,2,3,4,5`  

2. Add cells to register  
Input values: Register name (string), Cell name (string), Number of cell copies (integer), Initial strand coverings (strand-coordinate pairs, separated by semi-colons; can be left blank)  
Adds cells to the rightmost side of a register. If the register name doesn't exist, a new register will be created. All registers are considered to be in the same solution, and the simulator will apply the instructions to all of them.  
Multiple copies of the same cell with the same initial strand covers can be added.  
Example:  
`Register 1101`  
`Bit`  
`2`  
`One-first,0;One-second,2`  
This will add two Bit cells to Register 1101, where each one has a One-first strand at index 0, and a One-second strand at index 2. Coordinates are zero-indexed. Visually, the result would look like the following:  
`=>==>|=>==>`  
`01234|01234`

3. Add strand type  
Input values: Strand name (string), Domain names (comma-separated strings), Is complementary to the top strand (Y or N), Color hex code (6-digit hexadecimal number, #000000 by default)  
Example:  
`One-first`  
`0,1`  
`N`  
`#541074`  
The complementary property determines whether the strand binds to the bottom strand of the register (solid line strands in the SIMD||DNA paper) or to the top strand of the register (dashed line strands, only used in instructions.) The color code determines the strand's color in the SVG output file.

4. Add instruction  
Input values: Number of instruction strands (integer), Strand names (strings, entered one by one)  
Adds a single instruction. The simulator assumes an infinite number of copies of each instruction strand is added each time the simulation is run.  
Example:  
`2`  
`1A`  
`1B`  

5. Add cell-strand labels  
Input values: Cell name, coordinate-strand pairs, and label (comma-separated string)  
Adds a string label underneath a cell in the SVG output if its top strands follow a certain pattern. This allows the user to easily see the values encoded by a cell at any point in the simulation without needing to remember the exact strand encodings.  
Example:  
`Bit,0,One-first,2,One-second,One`

6. Run simulation  
Applies the current instruction sequence to all registers. An SVG file is automatically output for each register in the solution. The final result is output at the end after all the instructions.

7. Save data  
Input value: File name (string)  
Saves the current data in a JSON format.

8. Turn step-by-step simulation on  
Default value: Off  
When turned on, the user will need to press Enter in the terminal after each instruction. While an ASCII representation is printed in the terminal after each step, the SVG file is only produced after the last instruction.

9. Keep results after simulation  
Default value: Don't keep  
When set to keep, after running an instruction cycle on the registers, the output register contents will replace the initial contents. 
This is helpful if the user wants to simulate multiple instruction cycles on the register, such as repeatedly incrementing a binary number with the binary increment example. 
Otherwise, the initial register contents set by the user will never change when running instruction cycles.

10. Show unused instruction strands  
Default value: Don't show  
When set to show, the simulator will display instructions strands that would have otherwise attached to the register if there had been an open toehold and no competing strands had displaced them. 
Otherwise, only instruction strands that attached successfully are shown.

11. Compress SVG drawings  
Default value: Don't compress  
When set to compress, in the SVG file, only instruction numbers are shown on the left, the cell height is reduced, and the unused instruction strands occupy the same height as the successful instruction strands.
Unused instruction strands are displayed with a dashed line in this mode, while successful instructions are displayed with a solid line, regardless of whether they're complementary to the top or bottom strand. 
When turned off, the left label states "Instruction [number]", the cell height is increased, and the unused instruction strands are placed higher than the successful instructions. 
Top complementary strands are shown with a dashed line and bottom complementary strands are shown with a solid line, similar to the convention used in the SIMD||DNA paper.

12. Convert turingmachine.io Turing machine to SIMD register  
Input value: File name of turingmachine.io YAML file (string)  
Takes a Turing machine specification from http://turingmachine.io/ and creates the SIMD||DNA representation. This will erase all existing data to make way for the Turing machine information. 
The register's contents are initialized to the string input on the Turing machine's tape as specified in the YAML file. This will reject any Turing machine that uses symbols other than 0, 1, and blank. 
If the conversion is successful, the result can be immediately shown through the Run simulation option.

13. Draw inert instructions in SVG
Default value: Don't draw
By default, if an instruction does not affect the register at all, the simulator will not draw that instruction in the SVG file. If turned on, the full instruction set will be draw, where inert instructions have a red X notated underneath the instruction number.

14. Exit  
Fin
