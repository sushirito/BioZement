#include "global.h"

void check_file(const char *filename, std::ifstream &file)
{
//	file.open(filename, std::ifstream::in);
	if (!file.good()) {
		printf("ERROR: Unable to open file %s\n", filename);
		MPI_Finalize();
		exit(0);
	}
}

FILE* my_fopen(const char *filename, const char *mode)
/*
my_fopen :
	opens a file and checks if it is open correctly

  INPUT  : 1) filename : name of the file
		   2) mode     : the file mode (read write et c.)

  OUTPUT : (FILE*) pointer to the opened file
*/ 
{
	FILE *fptr; /* file pointer */

	if ( (fptr = fopen(filename, mode)) == NULL )
	{
		printf("ERROR: Unable to open file %s\n", filename);
		MPI_Finalize();
		exit(0);
    }	

	return fptr;
}

/****************************************************************** getword */
int getword(char *s, char *t)
{
  static char *sl;
  int i;

  if (s)
    sl = s;
  for (i=0; (t[i] = *sl) != ' ' && t[i] != '\0'; i++, sl++)
    ;
  if (t[i] == ' ') {
    t[i] = '\0';
    sl++;
  }
  return i;
}



/****************************************************************** isdbl */
int isdbl(char *s)
{
  if (*s == '-' || *s == '+')
    s++;
  while (*s >= '0' && *s<='9')
    s++;
  if (*s == '.')
    while (*++s >= '0' && *s<='9')
      ;
  if ((*s == 'e' || *s == 'E') && *++s != '\0')
    return isint(s);
  return (*s == '\0') ? 1 : 0;
}
/****************************************************************** atodbl */
//double atodbl(char *s)
//{
//  double d, ee, pow=1.0, des, sign = (*s == '-') ? -1.0 : 1.0;
//  int n, num_read = 0;
//
//  if (*s == '-' || *s == '+')
//    s++;
//  for (d=0.0; *s>='0' && *s<='9'; s++) {
//    d = 10.0*d + *s - '0';
//    num_read = 1;
//  }
//  des = 1.0;
//  if (*s == '.') {
//    for (*s++; *s>='0' && *s<='9'; s++, des *= 10.0)
//      d = 10.0*d + *s - '0';
//    num_read = 1;
//  }
//  if (*s == 'e' || *s == 'E') {
//    if (!num_read)
//      d = 1.0;
//    if ((n=atoint(++s)) < 0) {
//      ee = 0.1;
//      n  = -n;
//    } else
//      ee = 10.0;
//    for (; n > 0; n--)
//      pow *= ee;
//  }
//  return sign * d * pow / des;
//}
/****************************************************************** isint */
int isint(char *s)
{
  if (*s == '+' || *s == '-')
    s++;
  if (*s>'9' || *s++<'0')
    return 0;
  while (*s<='9' && *s>='0')
    s++;
  if (*s == '\0')
    return 1;
  return 0;
}
/****************************************************************** atoint */
//int atoint(char *s)
//{
//  int n, sign = (*s == '-') ? -1 : 1;
//
//  if (*s == '-' || *s == '+')
//    s++;
//  for (n=0; *s>='0' && *s<='9'; *s++)
//    n = 10*n + *s - '0';
//  return sign*n;
//}
/****************************************************************** INPUT FILE FUNCTIONS */
/****************************************************************** fgetline */
int fgetline(FILE *fp, char s[])
{
  int i, c;

  do {
    i = 0;
    while ((c = fgetc(fp)) != EOF && c != '\r' && c != '\n')
      s[i++] = c;
    if (c == '\r') {
      c = fgetc(fp);
      /*      if (c == '\n') {
	      printf("FEIL\n");
	      exit(0);
	      } */
    }
    s[i] = '\0';
    cleanline(s);
  } while ( c != EOF && s[0] == '\0');
  return c;
}

/****************************************************************** cleanline */
int cleanline(char s[])
{
  int i, j, in_comment = 0, in_space = 1;


  for (i=0, j=0; s[i] != '\0'; i++) {
	  if (s[i] == '/' && s[i+1] == '*') {
      i++;
      in_comment = 1;
    }
    else if (s[i] == '*' && s[i+1] == '/') {
      if (!in_comment) {
	printf("error: no comment to end.\n");
	return 0;
      } else {
	i++;
	in_comment = 0;
      }
    }
    else if (!in_comment && (s[i] == ' ' || s[i] == '\t')) {
      if (!in_space) {
	in_space = 1;
	s[j++] = ' ';
      }
    }
    else if (!in_comment) {
      s[j++] = s[i]; /*(s[i] >= 'A' && s[i] <= 'Z') ? (s[i]-'A'+'a') : s[i]; */
      in_space = 0;
    }
  }
  if (in_comment) {
    printf("error: line ended before end of comment.\n");
    return 0;
  }
  if (j>0 && s[j-1] == ' ')
    s[j-1] = '\0';
  else
    s[j] = '\0';

  if(s[0] == '#') s[0] = '\0';

  return 1;
}


