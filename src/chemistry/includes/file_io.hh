/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
#ifndef AMANZI_CHEMISTRY_FILEIO_HH_
#define AMANZI_CHEMISTRY_FILEIO_HH_

#define MAXCARDLENGTH 5 // add 1 to account of end of line \0
#define MAXWORDLENGTH 33
#define MAXSTRINGLENGTH 1025

#include <cstdlib>
#include <cstring>

#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

class FileIO {

 public:

  FileIO(char *filename);
  FileIO(string filename);
  virtual ~FileIO();

  int getLine();
  int getInputLine();
  int readDouble(double *d);
  double readDoubleFast();
  int readInt(int *i);
  int readWord(char *word);
  int readQuotedWords(char *words);
  static int removeQuotes(char *str);
  int comparesTo(char *str);
  int startsWith(char *str);
  int findStringInFile(char *card);
  static void checkDefaultMessage(char *word, int *ierr);
  static void checkErrorMessage(char *word1, char *word2, int ierr);
  static void checkLineErrorMessage(char *word, int ierr);
  static void toLower(char *word);
  static void toUpper(char *word);

  fstream file;
  stringstream *buffer;

};

#endif // AMANZI_CHEMISTRY_FILEIO_HH_
