//
// Created by user on 13.10.19.
//

#ifndef SYSTEM_SOLVER__LATEXWRITER_H_
#define SYSTEM_SOLVER__LATEXWRITER_H_

class LatexWriter {
 private:
  string filename_;
  string START_STR = "\\documentclass[a4paper,12pt,fleqn]{article}\n"
                     "\\usepackage{amsmath}\n"
                     "\\usepackage{systeme}\n"
                     "\\usepackage[russian]{babel}\n"
                     "\n"
                     "\\begin{document}\n"
                     "\n"
                     "\\author{It's me, Le Chat}\n"
                     "\\title{Hello Latex}\n"
                     "\\maketitle";
  string END_STR = "\\end{document}";
  ofstream stream;

 public:
  explicit LatexWriter(const string& filename="latex.tex") : filename_(filename), stream(filename_, std::ios::out) {
    StartFile();
  }

  /*void Testing() {
    Matrix<int> A = generateRandomSymmetricMatrix<int>(3);
    this->WriteMatrix(A);
    EndFile();
    this->CompileThisHell();
  }*/

//  ~LatexWriter() {
//    EndFile();
//  }

  void PrintText(const string& text) {
    stream << text << endl;
  }

  void PrintFormula(const string& text) {
    stream << "$" << endl;
    stream << text << endl;
    stream << "$" << endl;
  }

  void CompileThisHell() {
    stream.close();
    system("latexmk -pdf latex.tex\n");
    stream.open(filename_, std::ios::out);
  }

  template<class T>
  void PrintMatrix(const Matrix<T>& a, const string& message = "") {
    stream << message << endl;
    stream << "$\\begin{bmatrix}" << endl;
    for (auto& i : a) {
      for (auto& j : i) {
        stream << setw(6) << j << " & ";
      }
      stream << "\\\\" << endl;
    }
    stream << "\\end{bmatrix}$" << endl;
  }

  void StartFile() {
    PrintText(START_STR);
  }
  void EndFile() {
    PrintText(END_STR);
  }

};

#endif //SYSTEM_SOLVER__LATEXWRITER_H_
