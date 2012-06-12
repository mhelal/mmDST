#include "scores.h"


int PAM250 (char char1, char char2) {
	int score;
	switch (char1) {
		case 'C':
			switch (char2) {
				case 'C':
					score = 12;
					break;
				case 'S':
					score = 0;
					break;
				case 'T':
					score = -2;
					break;
				case 'P':
					score = -3;
					break;
				case 'A':
					score = -2;
					break;
				case 'G':
					score = -3;
					break;
				case 'N':
					score = -4;
					break;
				case 'D':
					score = -5;
					break;
				case 'E':
					score = -5;
					break;
				case 'Q':
					score = -5;
					break;
				case 'H':
					score = -3;
					break;
				case 'R':
					score = -4;
					break;
				case 'K':
					score = -5;
					break;
				case 'M':
					score = -5;
					break;
				case 'I':
					score = -2;
					break;
				case 'L':
					score = -6;
					break;
				case 'V':
					score = -2;
					break;
				case 'F':
					score = -4;
					break;
				case 'Y':
					score = 0;
					break;
				case 'W':
					score = -8;
					break;
				case 'B':
					score = -4;
					break;
				case 'Z':
					score = -5;
					break;
				case 'X':
					score = -3;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case 'S':
			switch (char2) {
				case 'C':
					score = 0;
					break;
				case 'S':
					score = 2;
					break;
				case 'T':
					score = 1;
					break;
				case 'P':
					score = 1;
					break;
				case 'A':
					score = 1;
					break;
				case 'G':
					score = 1;
					break;
				case 'N':
					score = 1;
					break;
				case 'D':
					score = 0;
					break;
				case 'E':
					score = 0;
					break;
				case 'Q':
					score = -1;
					break;
				case 'H':
					score = -1;
					break;
				case 'R':
					score = 0;
					break;
				case 'K':
					score = 0;
					break;
				case 'M':
					score = -2;
					break;
				case 'I':
					score = -1;
					break;
				case 'L':
					score = -3;
					break;
				case 'V':
					score = -1;
					break;
				case 'F':
					score = -3;
					break;
				case 'Y':
					score = -3;
					break;
				case 'W':
					score = -2;
					break;
				case 'B':
					score = 0;
					break;
				case 'Z':
					score = 0;
					break;
				case 'X':
					score = 0;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case 'T':
			switch (char2) {
				case 'C':
					score = -2;
					break;
				case 'S':
					score = 1;
					break;
				case 'T':
					score = 3;
					break;
				case 'P':
					score = 0;
					break;
				case 'A':
					score = 1;
					break;
				case 'G':
					score = 0;
					break;
				case 'N':
					score = 0;
					break;
				case 'D':
					score = 0;
					break;
				case 'E':
					score = 0;
					break;
				case 'Q':
					score = -1;
					break;
				case 'H':
					score = -1;
					break;
				case 'R':
					score = -1;
					break;
				case 'K':
					score = 0;
					break;
				case 'M':
					score = -1;
					break;
				case 'I':
					score = 0;
					break;
				case 'L':
					score = -2;
					break;
				case 'V':
					score = 0;
					break;
				case 'F':
					score = -3;
					break;
				case 'Y':
					score = -3;
					break;
				case 'W':
					score = -5;
					break;
				case 'B':
					score = 0;
					break;
				case 'Z':
					score = -1;
					break;
				case 'X':
					score = 0;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case 'P':
			switch (char2) {
				case 'C':
					score = -3;
					break;
				case 'S':
					score = 1;
					break;
				case 'T':
					score = 0;
					break;
				case 'P':
					score = 6;
					break;
				case 'A':
					score = 1;
					break;
				case 'G':
					score = -1;
					break;
				case 'N':
					score = -1;
					break;
				case 'D':
					score = -1;
					break;
				case 'E':
					score = -1;
					break;
				case 'Q':
					score = 0;
					break;
				case 'H':
					score = 0;
					break;
				case 'R':
					score = 0;
					break;
				case 'K':
					score = -1;
					break;
				case 'M':
					score = -2;
					break;
				case 'I':
					score = -2;
					break;
				case 'L':
					score = -3;
					break;
				case 'V':
					score = -1;
					break;
				case 'F':
					score = -5;
					break;
				case 'Y':
					score = -5;
					break;
				case 'W':
					score = -6;
					break;
				case 'B':
					score = -1;
					break;
				case 'Z':
					score = 0;
					break;
				case 'X':
					score = -1;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case 'A':
			switch (char2) {
				case 'C':
					score = -2;
					break;
				case 'S':
					score = 1;
					break;
				case 'T':
					score = 1;
					break;
				case 'P':
					score = 1;
					break;
				case 'A':
					score = 2;
					break;
				case 'G':
					score = 1;
					break;
				case 'N':
					score = 0;
					break;
				case 'D':
					score = 0;
					break;
				case 'E':
					score = 0;
					break;
				case 'Q':
					score = 0;
					break;
				case 'H':
					score = -1;
					break;
				case 'R':
					score = -2;
					break;
				case 'K':
					score = -1;
					break;
				case 'M':
					score = -1;
					break;
				case 'I':
					score = -1;
					break;
				case 'L':
					score = -2;
					break;
				case 'V':
					score = -0;
					break;
				case 'F':
					score = -4;
					break;
				case 'Y':
					score = -3;
					break;
				case 'W':
					score = -6;
					break;
				case 'B':
					score = 0;
					break;
				case 'Z':
					score = 0;
					break;
				case 'X':
					score = 0;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case 'G':
			switch (char2) {
				case 'C':
					score = -3;
					break;
				case 'S':
					score = 1;
					break;
				case 'T':
					score = 0;
					break;
				case 'P':
					score = -1;
					break;
				case 'A':
					score = 1;
					break;
				case 'G':
					score = 5;
					break;
				case 'N':
					score = 0;
					break;
				case 'D':
					score = 1;
					break;
				case 'E':
					score = 0;
					break;
				case 'Q':
					score = -1;
					break;
				case 'H':
					score = -2;
					break;
				case 'R':
					score = -3;
					break;
				case 'K':
					score = -2;
					break;
				case 'M':
					score = -3;
					break;
				case 'I':
					score = -3;
					break;
				case 'L':
					score = -4;
					break;
				case 'V':
					score = -1;
					break;
				case 'F':
					score = -5;
					break;
				case 'Y':
					score = -5;
					break;
				case 'W':
					score = -7;
					break;
				case 'B':
					score = 0;
					break;
				case 'Z':
					score = 0;
					break;
				case 'X':
					score = -1;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case 'N':
			switch (char2) {
				case 'C':
					score = -4;
					break;
				case 'S':
					score = 1;
					break;
				case 'T':
					score = 0;
					break;
				case 'P':
					score = -1;
					break;
				case 'A':
					score = 0;
					break;
				case 'G':
					score = 0;
					break;
				case 'N':
					score = 2;
					break;
				case 'D':
					score = 2;
					break;
				case 'E':
					score = 1;
					break;
				case 'Q':
					score = 1;
					break;
				case 'H':
					score = 2;
					break;
				case 'R':
					score = 0;
					break;
				case 'K':
					score = 1;
					break;
				case 'M':
					score = -2;
					break;
				case 'I':
					score = -2;
					break;
				case 'L':
					score = -3;
					break;
				case 'V':
					score = -2;
					break;
				case 'F':
					score = -4;
					break;
				case 'Y':
					score = -2;
					break;
				case 'W':
					score = -4;
					break;
				case 'B':
					score = 2;
					break;
				case 'Z':
					score = 2;
					break;
				case 'X':
					score = -1;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case 'D':
			switch (char2) {
				case 'C':
					score = -5;
					break;
				case 'S':
					score = 0;
					break;
				case 'T':
					score = 0;
					break;
				case 'P':
					score = -1;
					break;
				case 'A':
					score = 0;
					break;
				case 'G':
					score = 1;
					break;
				case 'N':
					score = 2;
					break;
				case 'D':
					score = 4;
					break;
				case 'E':
					score = 3;
					break;
				case 'Q':
					score = 2;
					break;
				case 'H':
					score = 1;
					break;
				case 'R':
					score = -1;
					break;
				case 'K':
					score = 0;
					break;
				case 'M':
					score = -3;
					break;
				case 'I':
					score = -2;
					break;
				case 'L':
					score = -4;
					break;
				case 'V':
					score = -2;
					break;
				case 'F':
					score = -6;
					break;
				case 'Y':
					score = -4;
					break;
				case 'W':
					score = -7;
					break;
				case 'B':
					score = 3;
					break;
				case 'Z':
					score = 3;
					break;
				case 'X':
					score = -1;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case 'E':
			switch (char2) {
				case 'C':
					score = -5;
					break;
				case 'S':
					score = 0;
					break;
				case 'T':
					score = 0;
					break;
				case 'P':
					score = -1;
					break;
				case 'A':
					score = 0;
					break;
				case 'G':
					score = 0;
					break;
				case 'N':
					score = 1;
					break;
				case 'D':
					score = 3;
					break;
				case 'E':
					score = 5;
					break;
				case 'Q':
					score = 2;
					break;
				case 'H':
					score = 1;
					break;
				case 'R':
					score = -1;
					break;
				case 'K':
					score = 0;
					break;
				case 'M':
					score = -2;
					break;
				case 'I':
					score = -2;
					break;
				case 'L':
					score = -3;
					break;
				case 'V':
					score = -2;
					break;
				case 'F':
					score = -5;
					break;
				case 'Y':
					score = -4;
					break;
				case 'W':
					score = -7;
					break;
				case 'B':
					score = 3;
					break;
				case 'Z':
					score = 3;
					break;
				case 'X':
					score = -1;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case 'Q':
			switch (char2) {
				case 'C':
					score = -5;
					break;
				case 'S':
					score = -1;
					break;
				case 'T':
					score = -1;
					break;
				case 'P':
					score = 0;
					break;
				case 'A':
					score = 0;
					break;
				case 'G':
					score = -1;
					break;
				case 'N':
					score = 1;
					break;
				case 'D':
					score = 2;
					break;
				case 'E':
					score = 2;
					break;
				case 'Q':
					score = 4;
					break;
				case 'H':
					score = 3;
					break;
				case 'R':
					score = 1;
					break;
				case 'K':
					score = 1;
					break;
				case 'M':
					score = -1;
					break;
				case 'I':
					score = -2;
					break;
				case 'L':
					score = -2;
					break;
				case 'V':
					score = -2;
					break;
				case 'F':
					score = -5;
					break;
				case 'Y':
					score = -4;
					break;
				case 'W':
					score = -5;
					break;
				case 'B':
					score = 1;
					break;
				case 'Z':
					score = 1;
					break;
				case 'X':
					score = -1;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case 'H':
			switch (char2) {
				case 'C':
					score = -3;
					break;
				case 'S':
					score = -1;
					break;
				case 'T':
					score = -1;
					break;
				case 'P':
					score = 0;
					break;
				case 'A':
					score = -1;
					break;
				case 'G':
					score = -2;
					break;
				case 'N':
					score = 2;
					break;
				case 'D':
					score = 1;
					break;
				case 'E':
					score = 1;
					break;
				case 'Q':
					score = 3;
					break;
				case 'H':
					score = 6;
					break;
				case 'R':
					score = 2;
					break;
				case 'K':
					score = 0;
					break;
				case 'M':
					score = -2;
					break;
				case 'I':
					score = -2;
					break;
				case 'L':
					score = -2;
					break;
				case 'V':
					score = -2;
					break;
				case 'F':
					score = -2;
					break;
				case 'Y':
					score = 0;
					break;
				case 'W':
					score = -3;
					break;
				case 'B':
					score = 1;
					break;
				case 'Z':
					score = 1;
					break;
				case 'X':
					score = -1;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case 'R':
			switch (char2) {
				case 'C':
					score = -4;
					break;
				case 'S':
					score = 0;
					break;
				case 'T':
					score = -1;
					break;
				case 'P':
					score = 0;
					break;
				case 'A':
					score = -2;
					break;
				case 'G':
					score = -3;
					break;
				case 'N':
					score = 0;
					break;
				case 'D':
					score = -1;
					break;
				case 'E':
					score = -1;
					break;
				case 'Q':
					score = 1;
					break;
				case 'H':
					score = 2;
					break;
				case 'R':
					score = 6;
					break;
				case 'K':
					score = 3;
					break;
				case 'M':
					score = 0;
					break;
				case 'I':
					score = -2;
					break;
				case 'L':
					score = -3;
					break;
				case 'V':
					score = -2;
					break;
				case 'F':
					score = -4;
					break;
				case 'Y':
					score = -4;
					break;
				case 'W':
					score = 2;
					break;
				case 'B':
					score = -1;
					break;
				case 'Z':
					score = -1;
					break;
				case 'X':
					score = -1;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case 'K':
			switch (char2) {
				case 'C':
					score = -5;
					break;
				case 'S':
					score = 0;
					break;
				case 'T':
					score = 0;
					break;
				case 'P':
					score = -1;
					break;
				case 'A':
					score = -1;
					break;
				case 'G':
					score = -2;
					break;
				case 'N':
					score = 1;
					break;
				case 'D':
					score = 0;
					break;
				case 'E':
					score = 0;
					break;
				case 'Q':
					score = 1;
					break;
				case 'H':
					score = 0;
					break;
				case 'R':
					score = 3;
					break;
				case 'K':
					score = 5;
					break;
				case 'M':
					score = 0;
					break;
				case 'I':
					score = -2;
					break;
				case 'L':
					score = -3;
					break;
				case 'V':
					score = -2;
					break;
				case 'F':
					score = -5;
					break;
				case 'Y':
					score = -4;
					break;
				case 'W':
					score = -3;
					break;
				case 'B':
					score = 1;
					break;
				case 'Z':
					score = 1;
					break;
				case 'X':
					score = -1;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case 'M':
			switch (char2) {
				case 'C':
					score = -5;
					break;
				case 'S':
					score = -2;
					break;
				case 'T':
					score = -1;
					break;
				case 'P':
					score = -2;
					break;
				case 'A':
					score = -1;
					break;
				case 'G':
					score = -3;
					break;
				case 'N':
					score = -2;
					break;
				case 'D':
					score = -3;
					break;
				case 'E':
					score = -2;
					break;
				case 'Q':
					score = -1;
					break;
				case 'H':
					score = -2;
					break;
				case 'R':
					score = 0;
					break;
				case 'K':
					score = 0;
					break;
				case 'M':
					score = 6;
					break;
				case 'I':
					score = 2;
					break;
				case 'L':
					score = 4;
					break;
				case 'V':
					score = 2;
					break;
				case 'F':
					score = 0;
					break;
				case 'Y':
					score = -2;
					break;
				case 'W':
					score = -4;
					break;
				case 'B':
					score = -2;
					break;
				case 'Z':
					score = -2;
					break;
				case 'X':
					score = -1;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case 'I':
			switch (char2) {
				case 'C':
					score = -2;
					break;
				case 'S':
					score = -1;
					break;
				case 'T':
					score = 0;
					break;
				case 'P':
					score = -2;
					break;
				case 'A':
					score = -1;
					break;
				case 'G':
					score = -3;
					break;
				case 'N':
					score = -2;
					break;
				case 'D':
					score = -2;
					break;
				case 'E':
					score = -2;
					break;
				case 'Q':
					score = -2;
					break;
				case 'H':
					score = -2;
					break;
				case 'R':
					score = -2;
					break;
				case 'K':
					score = -2;
					break;
				case 'M':
					score = 2;
					break;
				case 'I':
					score = 5;
					break;
				case 'L':
					score = 2;
					break;
				case 'V':
					score = 4;
					break;
				case 'F':
					score = 1;
					break;
				case 'Y':
					score = -1;
					break;
				case 'W':
					score = -5;
					break;
				case 'B':
					score = -2;
					break;
				case 'Z':
					score = -2;
					break;
				case 'X':
					score = -1;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case 'L':
			switch (char2) {
				case 'C':
					score = -6;
					break;
				case 'S':
					score = -3;
					break;
				case 'T':
					score = -2;
					break;
				case 'P':
					score = -3;
					break;
				case 'A':
					score = -2;
					break;
				case 'G':
					score = -4;
					break;
				case 'N':
					score = -3;
					break;
				case 'D':
					score = -4;
					break;
				case 'E':
					score = -3;
					break;
				case 'Q':
					score = -2;
					break;
				case 'H':
					score = -2;
					break;
				case 'R':
					score = -3;
					break;
				case 'K':
					score = -3;
					break;
				case 'M':
					score = 4;
					break;
				case 'I':
					score = 2;
					break;
				case 'L':
					score = 6;
					break;
				case 'V':
					score = 2;
					break;
				case 'F':
					score = 2;
					break;
				case 'Y':
					score = -1;
					break;
				case 'W':
					score = -2;
					break;
				case 'B':
					score = -3;
					break;
				case 'Z':
					score = -3;
					break;
				case 'X':
					score = -1;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case 'V':
			switch (char2) {
				case 'C':
					score = -2;
					break;
				case 'S':
					score = -1;
					break;
				case 'T':
					score = 0;
					break;
				case 'P':
					score = -1;
					break;
				case 'A':
					score = 0;
					break;
				case 'G':
					score = -1;
					break;
				case 'N':
					score = -2;
					break;
				case 'D':
					score = -2;
					break;
				case 'E':
					score = -2;
					break;
				case 'Q':
					score = -2;
					break;
				case 'H':
					score = -2;
					break;
				case 'R':
					score = -2;
					break;
				case 'K':
					score = -2;
					break;
				case 'M':
					score = 2;
					break;
				case 'I':
					score = 4;
					break;
				case 'L':
					score = 2;
					break;
				case 'V':
					score = 4;
					break;
				case 'F':
					score = -1;
					break;
				case 'Y':
					score = -2;
					break;
				case 'W':
					score = -6;
					break;
				case 'B':
					score = -2;
					break;
				case 'Z':
					score = -2;
					break;
				case 'X':
					score = -1;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case 'F':
			switch (char2) {
				case 'C':
					score = -4;
					break;
				case 'S':
					score = -3;
					break;
				case 'T':
					score = -3;
					break;
				case 'P':
					score = -5;
					break;
				case 'A':
					score = -4;
					break;
				case 'G':
					score = -5;
					break;
				case 'N':
					score = -4;
					break;
				case 'D':
					score = -6;
					break;
				case 'E':
					score = -5;
					break;
				case 'Q':
					score = -5;
					break;
				case 'H':
					score = -2;
					break;
				case 'R':
					score = -4;
					break;
				case 'K':
					score = -5;
					break;
				case 'M':
					score = 0;
					break;
				case 'I':
					score = 1;
					break;
				case 'L':
					score = 2;
					break;
				case 'V':
					score = -1;
					break;
				case 'F':
					score = 9;
					break;
				case 'Y':
					score = 7;
					break;
				case 'W':
					score = 0;
					break;
				case 'B':
					score = -4;
					break;
				case 'Z':
					score = -5;
					break;
				case 'X':
					score = -2;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case 'Y':
			switch (char2) {
				case 'C':
					score = 0;
					break;
				case 'S':
					score = -3;
					break;
				case 'T':
					score = -3;
					break;
				case 'P':
					score = -5;
					break;
				case 'A':
					score = -3;
					break;
				case 'G':
					score = -5;
					break;
				case 'N':
					score = -2;
					break;
				case 'D':
					score = -4;
					break;
				case 'E':
					score = -4;
					break;
				case 'Q':
					score = -4;
					break;
				case 'H':
					score = 0;
					break;
				case 'R':
					score = -4;
					break;
				case 'K':
					score = -4;
					break;
				case 'M':
					score = -2;
					break;
				case 'I':
					score = -1;
					break;
				case 'L':
					score = -1;
					break;
				case 'V':
					score = -2;
					break;
				case 'F':
					score = 7;
					break;
				case 'Y':
					score = 10;
					break;
				case 'W':
					score = 0;
					break;
				case 'B':
					score = -3;
					break;
				case 'Z':
					score = -4;
					break;
				case 'X':
					score = -2;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case 'W':
			switch (char2) {
				case 'C':
					score = -8;
					break;
				case 'S':
					score = -2;
					break;
				case 'T':
					score = -5;
					break;
				case 'P':
					score = -6;
					break;
				case 'A':
					score = -6;
					break;
				case 'G':
					score = -7;
					break;
				case 'N':
					score = -4;
					break;
				case 'D':
					score = -7;
					break;
				case 'E':
					score = -7;
					break;
				case 'Q':
					score = -5;
					break;
				case 'H':
					score = -3;
					break;
				case 'R':
					score = 2;
					break;
				case 'K':
					score = -3;
					break;
				case 'M':
					score = -4;
					break;
				case 'I':
					score = -5;
					break;
				case 'L':
					score = -2;
					break;
				case 'V':
					score = -6;
					break;
				case 'F':
					score = 0;
					break;
				case 'Y':
					score = 0;
					break;
				case 'W':
					score = 17;
					break;
				case 'B':
					score = -5;
					break;
				case 'Z':
					score = -6;
					break;
				case 'X':
					score = -4;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case 'B':
			switch (char2) {
				case 'C':
					score = -4;
					break;
				case 'S':
					score = 0;
					break;
				case 'T':
					score = 0;
					break;
				case 'P':
					score = -1;
					break;
				case 'A':
					score = 0;
					break;
				case 'G':
					score = 0;
					break;
				case 'N':
					score = 2;
					break;
				case 'D':
					score = 3;
					break;
				case 'E':
					score = 3;
					break;
				case 'Q':
					score = 1;
					break;
				case 'H':
					score = 1;
					break;
				case 'R':
					score = -1;
					break;
				case 'K':
					score = 1;
					break;
				case 'M':
					score = -2;
					break;
				case 'I':
					score = -2;
					break;
				case 'L':
					score = -3;
					break;
				case 'V':
					score = -2;
					break;
				case 'F':
					score = -4;
					break;
				case 'Y':
					score = -3;
					break;
				case 'W':
					score = -5;
					break;
				case 'B':
					score = 3;
					break;
				case 'Z':
					score = 2;
					break;
				case 'X':
					score = -1;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case 'Z':
			switch (char2) {
				case 'C':
					score = -5;
					break;
				case 'S':
					score = 0;
					break;
				case 'T':
					score = -1;
					break;
				case 'P':
					score = 0;
					break;
				case 'A':
					score = 0;
					break;
				case 'G':
					score = 0;
					break;
				case 'N':
					score = 2;
					break;
				case 'D':
					score = 3;
					break;
				case 'E':
					score = 3;
					break;
				case 'Q':
					score = 1;
					break;
				case 'H':
					score = 1;
					break;
				case 'R':
					score = -1;
					break;
				case 'K':
					score = 1;
					break;
				case 'M':
					score = -2;
					break;
				case 'I':
					score = -2;
					break;
				case 'L':
					score = -3;
					break;
				case 'V':
					score = -2;
					break;
				case 'F':
					score = -5;
					break;
				case 'Y':
					score = -4;
					break;
				case 'W':
					score = -6;
					break;
				case 'B':
					score = 2;
					break;
				case 'Z':
					score = 3;
					break;
				case 'X':
					score = -1;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case 'X':
			switch (char2) {
				case 'C':
					score = -3;
					break;
				case 'S':
					score = 0;
					break;
				case 'T':
					score = 0;
					break;
				case 'P':
					score = -1;
					break;
				case 'A':
					score = 0;
					break;
				case 'G':
					score = -1;
					break;
				case 'N':
					score = -1;
					break;
				case 'D':
					score = -1;
					break;
				case 'E':
					score = -1;
					break;
				case 'Q':
					score = -1;
					break;
				case 'H':
					score = -1;
					break;
				case 'R':
					score = -1;
					break;
				case 'K':
					score = -1;
					break;
				case 'M':
					score = -1;
					break;
				case 'I':
					score = -1;
					break;
				case 'L':
					score = -1;
					break;
				case 'V':
					score = -1;
					break;
				case 'F':
					score = -2;
					break;
				case 'Y':
					score = -2;
					break;
				case 'W':
					score = -4;
					break;
				case 'B':
					score = -1;
					break;
				case 'Z':
					score = -1;
					break;
				case 'X':
					score = -1;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
		case '-':
			switch (char2) {
				case 'C':
					score = -8;
					break;
				case 'S':
					score = -8;
					break;
				case 'T':
					score = -8;
					break;
				case 'P':
					score = -8;
					break;
				case 'A':
					score = -8;
					break;
				case 'G':
					score = -8;
					break;
				case 'N':
					score = -8;
					break;
				case 'D':
					score = -8;
					break;
				case 'E':
					score = -8;
					break;
				case 'Q':
					score = -8;
					break;
				case 'H':
					score = -8;
					break;
				case 'R':
					score = -8;
					break;
				case 'K':
					score = -8;
					break;
				case 'M':
					score = -8;
					break;
				case 'I':
					score = -8;
					break;
				case 'L':
					score = -8;
					break;
				case 'V':
					score = -8;
					break;
				case 'F':
					score = -8;
					break;
				case 'Y':
					score = -8;
					break;
				case 'W':
					score = -8;
					break;
				case 'B':
					score = -8;
					break;
				case 'Z':
					score = -8;
					break;
				case 'X':
					score = -8;
					break;
				case '-':
					score = -8;
					break;
			}
			break;
	}
	return score;
}

int BLOSUM (char char1, char char2) {
	int score;
	switch (char1) {
		case 'A':
			switch (char2) {
				case 'A':
					score = 1;
					break;
				case 'C':
					score = -1;
					break;
				case 'G':
					score = -1;
					break;
				case 'T':
					score = -1;
					break;
				case 'U':
					score = -1;
					break;
				case 'X':
					score = -1;
					break;
				case '-':
					score = -2;
					break;
			}
			break;
		case 'C':
			switch (char2) {
				case 'A':
					score = -1;
					break;
				case 'C':
					score = 1;
					break;
				case 'G':
					score = -1;
					break;
				case 'T':
					score = -1;
					break;
				case 'U':
					score = -1;
					break;
				case 'X':
					score = -1;
					break;
				case '-':
					score = -2;
					break;
			}
			break;
		case 'G':
			switch (char2) {
				case 'A':
					score = -1;
					break;
				case 'C':
					score = -1;
					break;
				case 'G':
					score = 1;
					break;
				case 'T':
					score = -1;
					break;
				case 'U':
					score = -1;
					break;
				case 'X':
					score = -1;
					break;
				case '-':
					score = -2;
					break;
			}
			break;
		case 'T':
			switch (char2) {
				case 'A':
					score = -1;
					break;
				case 'C':
					score = -1;
					break;
				case 'G':
					score = -1;
					break;
				case 'T':
					score = 1;
					break;
				case 'U':
					score = -1;
					break;
				case 'X':
					score = -1;
					break;
				case '-':
					score = -2;
					break;
			}
			break;
		case 'U':
			switch (char2) {
				case 'A':
					score = -1;
					break;
				case 'C':
					score = -1;
					break;
				case 'G':
					score = -1;
					break;
				case 'T':
					score = -1;
					break;
				case 'U':
					score = 1;
					break;
				case 'X':
					score = -1;
					break;
				case '-':
					score = -2;
					break;
			}
			break;
		case 'X':
			switch (char2) {
				case 'A':
					score = -1;
					break;
				case 'C':
					score = -1;
					break;
				case 'G':
					score = -1;
					break;
				case 'T':
					score = -1;
					break;
				case 'U':
					score = -1;
					break;
				case 'X':
					score = 1;
					break;
				case '-':
					score = -2;
					break;
			}
			break;
	}
	return score;
}
