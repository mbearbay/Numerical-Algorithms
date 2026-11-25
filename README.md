# Numerical Algorithms

[English](#english) | [Dansk](#dansk)

## English

Exam project from the course *Numerical Algorithms* at Copenhagen Business School, written in C++ without external libraries. The program fits regression splines to data using least squares.

### What the program does

1. Reads x and y values from files.
2. Asks the user for the number of intervals and the order of the spline.
3. Builds the knot vector.
4. Evaluates the B-splines with **de Boor's algorithm** and builds the system matrix.
5. Solves the least squares problem with **Householder QR factorisation** and back substitution.
6. Prints the spline coefficients.

All vector and matrix operations (dot product, norm, QR factorisation) are implemented from scratch.

### Files

| File | Content |
| --- | --- |
| `Code NA eksamen.cpp` | The main program |
| `Code for external tests` | Extra test functions, e.g. Hilbert matrices and normal equations |
| `Routediagram for code.png` | Diagram of how the functions are connected |
| `NA_Eksamensbeskrivelse.pdf` | The exam assignment |

### How to run

```bash
g++ -std=c++17 -O2 "Code NA eksamen.cpp" -o spline
./spline
```

The comments in the code and the program's text are in Danish.

**Tools:** C++ (standard library only)

---

## Dansk

Eksamensprojekt fra faget *Numeriske algoritmer* på Copenhagen Business School, skrevet i C++ uden eksterne biblioteker. Programmet tilpasser regressionssplines til data med mindste kvadraters metode.

### Hvad programmet gør

1. Indlæser x- og y-værdier fra filer.
2. Spørger brugeren om antal intervaller og splinens orden.
3. Danner knudepunkterne.
4. Evaluerer B-splines med **De Boors algoritme** og danner totalmatricen.
5. Løser mindste kvadraters problem med **Householders QR-faktorisering** og baglæns substitution.
6. Udskriver splinens koefficienter.

Alle vektor- og matrixoperationer (prikprodukt, norm, QR-faktorisering) er implementeret fra bunden.

### Filer

| Fil | Indhold |
| --- | --- |
| `Code NA eksamen.cpp` | Hovedprogrammet |
| `Code for external tests` | Ekstra testfunktioner, fx Hilbert-matricer og normalligninger |
| `Routediagram for code.png` | Diagram over hvordan funktionerne hænger sammen |
| `NA_Eksamensbeskrivelse.pdf` | Eksamensopgaven |

### Sådan køres koden

```bash
g++ -std=c++17 -O2 "Code NA eksamen.cpp" -o spline
./spline
```

**Værktøjer:** C++ (kun standardbiblioteket)
