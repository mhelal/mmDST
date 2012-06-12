#ifndef SCORESH_
#define SCORESH_

#define GAPCHAR '-' 
int gapOpenning;
int gapExtension;
int matchscore;
int mismatchscore;
char * m_csScoreMatrixName;
int m_csScoreMatrixCode;

int subScore (char char1, char char2, int stype);
int gapScore(int stype);
int GetProtein(char ch);
int getDNA (char dna);
char encodeDNA (int dna);
int getScore(char ch1,char ch2);
void SetScoreMatrix(char * str, int Matrix);
enum PROTEINS{
    A = 1,		// Alanine		(Ala)		Aliphatic side chains (these are hydrophobic - they like to be away from water)
    R = 2,		// Arginine		(Arg)		Very hydrophilic  
    N = 3,		// Asparagine		(Asn)		Acidic  
    D = 4,		// Aspartate		(Asp)		Acidic
    C = 5,		// Cysteine		(Cys)		Sulfur containing amino acids - can form disulfide bonds with each other
    Q = 6,		// Glutamine		(Gln)		Acidic 
    E = 7,		// Glutamate		(Glu)		Acidic
    G = 8,		// Glycine		(Gly)		Aliphatic side chains (these are hydrophobic - they like to be away from water)
    H = 9,		// Histidine		(His)		Very hydrophilic  
    I =10,		// Isoleucine		(Ile)		Aliphatic side chains (these are hydrophobic - they like to be away from water)  
    L =11,		// Leucine		(Leu)		Aliphatic side chains (these are hydrophobic - they like to be away from water)  
    K =12,		// Lysine		(Lys)		Very hydrophilic
    M =13,		// Methionine		(Met)		Sulfur containing amino acids - can form disulfide bonds with each other  
    F =14,		// Phenlylalanine	(Phe)		Aromatic side chains - more reactive than other amino acids
    P =15,		// Proline		(Pro)		Hydrophilic - in a category by itself, shape results in bending of the amino acid chain
    S =16,		// Serine		(Ser)		Hydroxyl group, somewhat hydrophilic
    T =17,		// Threonine		(Thr)		Hydroxyl group, somewhat hydrophilic
    W =18,		// Tryptophan		(Trp)		Aromatic side chains - more reactive than other amino acids  
    Y =19,		// Tyrosine		(Tyr)		Aromatic side chains - more reactive than other amino acids    
    V =20,		// Valine		(Val)		Aliphatic side chains (these are hydrophobic - they like to be away from water)  
    B =21,  
    Z =22,  
    X =23,  
    GAP=24
};


short * m_pScoreMatrix[24];
/*
static short m_BLOSUM30[24][24];
static short m_BLOSUM35[24][24];
static short m_BLOSUM40[24][24];
static short m_BLOSUM45[24][24];
static short m_BLOSUM50[24][24];
static short m_BLOSUM55[24][24];
static short m_BLOSUM60[24][24];
static short m_BLOSUM62[24][24];
static short m_BLOSUM65[24][24];
static short m_BLOSUM70[24][24];
static short m_BLOSUM75[24][24];
static short m_BLOSUM80[24][24];
static short m_BLOSUM85[24][24];
static short m_BLOSUM90[24][24];
static short m_BLOSUM100[24][24];
static short m_BLOSUMN[24][24];


static short m_PAM10[24][24];
static short m_PAM20[24][24];
static short m_PAM30[24][24];
static short m_PAM40[24][24];
static short m_PAM50[24][24];
static short m_PAM60[24][24];
static short m_PAM70[24][24];
static short m_PAM80[24][24];
static short m_PAM90[24][24];

static short m_PAM100[24][24];
static short m_PAM110[24][24];
static short m_PAM120[24][24];
static short m_PAM130[24][24];
static short m_PAM140[24][24];
static short m_PAM150[24][24];
static short m_PAM160[24][24];
static short m_PAM170[24][24];
static short m_PAM180[24][24];
static short m_PAM190[24][24];

static short m_PAM200[24][24];
static short m_PAM210[24][24];
static short m_PAM220[24][24];
static short m_PAM230[24][24];
static short m_PAM240[24][24];
static short m_PAM250[24][24];
static short m_PAM260[24][24];
static short m_PAM270[24][24];
static short m_PAM280[24][24];
static short m_PAM290[24][24];

static short m_PAM300[24][24];
static short m_PAM310[24][24];
static short m_PAM320[24][24];
static short m_PAM330[24][24];
static short m_PAM340[24][24];
static short m_PAM350[24][24];
static short m_PAM360[24][24];
static short m_PAM370[24][24];
static short m_PAM380[24][24];
static short m_PAM390[24][24];

static short m_PAM400[24][24];
static short m_PAM410[24][24];
static short m_PAM420[24][24];
static short m_PAM430[24][24];
static short m_PAM440[24][24];
static short m_PAM450[24][24];
static short m_PAM460[24][24];
static short m_PAM470[24][24];
static short m_PAM480[24][24];
static short m_PAM490[24][24];
static short m_PAM500[24][24];

static short m_MATCH[24][24];
*/
#endif
