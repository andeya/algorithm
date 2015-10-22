/*
Throwdice rolls virtual dice and prints the result.

It reads dice-descriptions from command-line arguments and rolls
each of them, printing the results to stdout.

For example:
    throwdice 3d6+1
could output:
    14 (3, 6, 4)
Where 14 is the total, and 3, 6, and 4 are the individual dice rolls.

The '-d6' flag can be passed to use the D6 system when rolling.

For example:
	throwdice -d6 4
could output:
	17 (3, 4, 2, 6, 2)
Where the 4th die, coming up "6", triggered an additional dice roll.

When invoked with no command-line arguments, it reads from
standart input instead, one description per line.  In this case,
it will stop reading when receiving EOF or the words "exit" or
"quit" (case insensitive).
*/
package documentation
