#include "scores.h"

typedef struct RBAlign_rec { /* ScoringData */
    MOATypeDimn seqNum;
    MOATypeShape * seqLen;
    char * * sequences;
    /*CArray	<CPole,CPole&>	m_ArrPoles,m_ArrPrimaryPoles,m_ArrStickyPoles;*/
    MOATypeShape * * * Poles;
    MOATypeShape * * * PrimaryPoles;
    MOATypeShape * * * StickyPoles;
    MOATypeShape PrimaryPolesNum;
    int FixType; /*Primary Data Fix type*/
    int m_NumOfRuns;
    int	m_SticknessPercentage;
    int	m_NumberOfPoles;
    int	m_NumberOfTries;
    int	m_TryCounter;
    int	m_RoundCounter;
    int	m_PPGeneration;
    int	m_PPSelection;
    int	m_NP;
    int	m_PPWideness;
    int	m_PPSimilarity;
    int	m_RunningMode;
    int	m_TryMode;    
} RBAlignData;

void FindInitialPrimaryPoles(RBAlignData RBData) {
    MOATypeDimn k, l;
    int i;
    MOATypeShape j, AC1,AC2;
    int	*pCounterLimits = NULL,*pCounters = NULL;
    int CorrFlag,SortFlag;
    CString str;
    MOATypeShape * * PrimaryPoles = NULL, * APrimaryPole = NULL;
    MOATypeShape * * * pAminoAcidIndeces = NULL;
    MOATypeShape * * pAminoAcidIndecesNum;
    MOATypeShape * * ArrDummyPoles = NULL;
    MOATypeShape ArrDummyPolesNum;
    int XX,MZ;
    double DD;
    int	MaxSize,MaxIndex,YY;
    int	MaxSeqLen;
    int	PPLocalWideness;
    int PPLocalSimilarity;

    PPLocalWideness=m_PPWideness+RandomInt(0,1);
    PPLocalSimilarity=m_PPSimilarity+RandomInt(0,1);

    pCounterLimits=malloc(RBData.seqNum * sizeof *int);
    pCounters=malloc(RBData.seqNum * sizeof *int);
    APrimaryPole=malloc(RBData.seqNum * sizeof *APrimaryPole);
    PrimaryPoles=malloc(RBData.seqNum * sizeof *PrimaryPoles);
    PrimaryPoles[0]=malloc(RBData.seqNum * sizeof *PrimaryPoles[0]);
    PrimaryPoles[1]=malloc(RBData.seqNum * sizeof *PrimaryPoles[1]);
    PrimaryPoles[2]=malloc(RBData.seqNum * sizeof *PrimaryPoles[2]);
    PrimaryPoles[3]=malloc(RBData.seqNum * sizeof *PrimaryPoles[3]);
    PrimaryPoles[4]=malloc(RBData.seqNum * sizeof *PrimaryPoles[4]);

    if (PrimaryPoles != NULL) {
        for (i=0;i<PrimaryPolesNum;i++)
            free(PrimaryPoles[i]);
        free (PrimaryPoles);
        PrimaryPolesNum = 0;
    }
    PrimaryPoles = NULL;

    pAminoAcidIndeces=malloc(RBData.seqNum * sizeof *pAminoAcidIndeces);
    pAminoAcidIndecesNum=malloc(RBData.seqNum * sizeof *pAminoAcidIndecesNum);
    for (i=0;i<RBData.seqNum;i++) {
        pAminoAcidIndeces[i]=malloc(20 * sizeof *pAminoAcidIndeces[i]);
        pAminoAcidIndecesNum[i] = malloc(20 * sizeof *pAminoAcidIndeces[i]);
        for (AC2=0;AC2<20;AC2++) 
            pAminoAcidIndecesNum[i][AC2] = 0;
    }

//Find Characters for each Array
    for (k=0;k<RBData.seqNum;k++) {
        for (j=0;j<RBData.seqLen[k];j++) {
            AC1=GetProtein(RBData.sequences[k][j])-1;
            for (AC2=0;AC2<20;AC2++) {
                if (m_pScoreMatrix[AC1][AC2]>=PPLocalSimilarity) {
                    pAminoAcidIndeces[k][AC2] = (MOATypeShape *) realloc (pAminoAcidIndeces[k][AC2], (pAminoAcidIndecesNum[k][AC2]+1) * sizeof *pAminoAcidIndeces[k][AC2]);
                    pAminoAcidIndeces[k][AC2][pAminoAcidIndecesNum[k][AC2]] = j;                    
                    pAminoAcidIndecesNum[k][AC2] ++;
                } //endif
            } //endfor
        } //endfor
    } //endfor


//Handle Empty Arrays
    if (m_PPGeneration==1 || m_PPGeneration==3)		{
        //Lenght Equalization - Find Less Probable Matches for Empty Arrays
        int Max,Min;
        for (i=0;i<20;i++) {
            Max=Min=int(pAminoAcidIndecesNum[0][i]);
            for (k=1;k<RBData.seqNum;k++) {
                    if (pAminoAcidIndecesNum[k][i]>Max)	Max=pAminoAcidIndecesNum[k][i];
                    if (pAminoAcidIndecesNum[k][i]<Min)	Min=pAminoAcidIndecesNum[k][i];
            } //endfor

            for (k=0;k<RBData.seqNum;k++) {
                if (DummyPoles != NULL) {
                    for (i=0;i<DummyPolesNum;i++)
                        free(DummyPoles[i]);
                    free (DummyPoles);
                    DummyPolesNum = 0;
                }
                DummyPoles = NULL;

                XX=PPLocalSimilarity;
                while (DummyPolesNum<Max && XX>-3) {
                    XX--;
                    for (j=0;j<20;j++) {
                        if (m_pScoreMatrix[i][j]==XX) {
                            DummyPoles = realloc (DummyPoles, (DummyPolesNum + 1) * sizeof *DummyPoles);
                            DummyPoles[DummyPolesNum] = malloc (BData.seqNum * sizeof * DummyPoles[DummyPolesNum]);
                            for (l=0;l<BData.seqNum;l++)
                                DummyPoles[DummyPolesNum][l] = pAminoAcidIndeces[k][j][l];
                            DummyPolesNum ++;
                        } //endif
                    } //endfor
                } //endwhile

                Min=Max*8/10;
                while (pAminoAcidIndecesNum[k][i]<Min && DummyPolesNum!=0) {
                        j=RandomInt(0,int(DummyPolesNum)-1);
                        pAminoAcidIndeces[k][i] = realloc (pAminoAcidIndeces[k][i], (pAminoAcidIndecesNum[k][i]+1) * sizeof *pAminoAcidIndeces[k][i]);                               
                        for (l=0;l<BData.seqNum;l++)
                            pAminoAcidIndeces[k][i][pAminoAcidIndecesNum[k][i]] = DummyPoles[j][l];
                        pAminoAcidIndecesNum[k][i] ++;
                        ArrDummyPoles.RemoveAt(j);
                } //endwhile
            } //endfor
        } //endfor
    } //endif

    // Delete Empty Rows
    for (i=0;i<20;i++) {
            CorrFlag=FALSE;
            for (k=0;k<RBData.seqNum;k++) {
                    if (pAminoAcidIndeces[k][i].GetSize()==0) {
                            CorrFlag=TRUE;
                            break;
                    } //endif
            } //endfor

            if (CorrFlag) {
                    for (k=0;k<RBData.seqNum;k++) {
                            pAminoAcidIndeces[k][i].RemoveAll();
                    } //endif
            } //endif
    } //endfor

//Sort Arrays
    for (k=0;k<RBData.seqNum;k++) {
            for (i=0;i<20;i++) {
                    SortFlag=FALSE;
                    while (!SortFlag) {
                            SortFlag=TRUE;
                            for (j=0;j<pAminoAcidIndeces[k][i].GetSize()-1;j++) {
                                    if (pAminoAcidIndeces[k][i][j]>pAminoAcidIndeces[k][i][j+1]) {
                                            XX=pAminoAcidIndeces[k][i][j];
                                            pAminoAcidIndeces[k][i][j]=pAminoAcidIndeces[k][i][j+1];
                                            pAminoAcidIndeces[k][i][j+1]=XX;

                                            SortFlag=FALSE;
                                    } //endif
                            } //endfor
                    } //endwhile
            } //endfor
    } //endfor


//	m_pRBSD->m_ExtraInfoList.AddString(_T("PP Generating"));

    if (m_pRBSD!=NULL) {
            m_pRBSD->m_InitializationProgreass.SetRange32(0,19);
            m_pRBSD->m_InitializationProgreass.SetPos(0);
    } //endif


//Generate Primary Poles by Index
    if (m_PPGeneration==0 || m_PPGeneration==1) {
            for (i=0;i<20;i++) { //All Amino Acids
                    HandleWindowsMessages();

                    if (m_pRBSD!=NULL)
                            m_pRBSD->m_InitializationProgreass.SetPos(i);

                    for (k=0;k<RBData.seqNum;k++) {
                            pCounterLimits[k]=int(pAminoAcidIndeces[k][i].GetSize());
                    } //endfor

                    MaxSize=pCounterLimits[0];
                    for (k=1;k<RBData.seqNum;k++) {
                            if (MaxSize<pCounterLimits[k])
                                    MaxSize=pCounterLimits[k];
                    } //endfor

                    for (j=0;j<MaxSize;j++) {
                            for (k=0;k<RBData.seqNum;k++) {
                                    MZ=MaxSize;
                                    if (MaxSize==1)	MZ++;

                                    XX=int(double(j*(pCounterLimits[k]-1))/(MZ-1)+0.5);
                                    PrimaryPoles[0][k]=pAminoAcidIndeces[k][i][XX];

                                    if (XX!=0)
                                            PrimaryPoles[1][k]=pAminoAcidIndeces[k][i][XX-1];
                                    else
                                            PrimaryPoles[1][k]=pAminoAcidIndeces[k][i][XX];

                                    if (XX!=pCounterLimits[k]-1)
                                            PrimaryPoles[2][k]=pAminoAcidIndeces[k][i][XX+1];
                                    else
                                            PrimaryPoles[2][k]=pAminoAcidIndeces[k][i][XX];

                                    if (XX!=0 && XX!=1)
                                            PrimaryPoles[3][k]=pAminoAcidIndeces[k][i][XX-2];
                                    else
                                            PrimaryPoles[3][k]=pAminoAcidIndeces[k][i][XX];

                                    if (XX!=pCounterLimits[k]-1 && XX!=pCounterLimits[k]-2)
                                            PrimaryPoles[4][k]=pAminoAcidIndeces[k][i][XX+2];
                                    else
                                            PrimaryPoles[4][k]=pAminoAcidIndeces[k][i][XX];
                            } //endfor

                            m_ArrPrimaryPoles.Add(PrimaryPoles[0]);

                            if (PPLocalWideness==1 || PPLocalWideness==2) {
                                    m_ArrPrimaryPoles.Add(PrimaryPoles[1]);
                                    m_ArrPrimaryPoles.Add(PrimaryPoles[2]);
                            } //endif

                            if (PPLocalWideness==2) {
                                    m_ArrPrimaryPoles.Add(PrimaryPoles[3]);
                                    m_ArrPrimaryPoles.Add(PrimaryPoles[4]);
                            } //endif

                    } //endfor
            } //endfor
    } //endif Index Generation


//Generate Primary Poles by Location
    if (m_PPGeneration==2 ||m_PPGeneration==3) {
            MaxSeqLen=m_InputSequences[0].m_Seq.GetLength();
            for (k=1;k<RBData.seqNum;k++) {
                    if (MaxSeqLen<m_InputSequences[k].m_Seq.GetLength())
                            MaxSeqLen=m_InputSequences[k].m_Seq.GetLength();
            } //endfor

            for (i=0;i<20;i++) { //All Amino Acids
                    HandleWindowsMessages();

                    for (k=0;k<RBData.seqNum;k++) {
                            pCounterLimits[k]=int(pAminoAcidIndeces[k][i].GetSize());
                    } //endfor

                    MaxSize=pCounterLimits[0];
                    MaxIndex=0;
                    for (k=1;k<RBData.seqNum;k++) {
                            if (MaxSize<pCounterLimits[k]) {
                                    MaxSize=pCounterLimits[k];
                                    MaxIndex=k;
                            } //endif
                    } //endfor

                    for (j=0;j<MaxSize;j++) {
                            SortFlag=TRUE;
                            for (k=0;k<RBData.seqNum;k++) {
                                    MZ=MaxSize;
                                    if (MaxSize==1)	MZ++;

                                    XX=pAminoAcidIndeces[MaxIndex][i][j]*m_InputSequences[k].m_Seq.GetLength()/m_InputSequences[MaxIndex].m_Seq.GetLength();
                                    YY=FindNearestPointInArray(XX,&pAminoAcidIndeces[k][i]);

                                    DD=double(YY-XX)/m_InputSequences[k].m_Seq.GetLength();
                                    if (abs(DD)<(PPLocalWideness*2+1)*0.05) {
                                            PrimaryPoles[0][k]=YY;
                                    } else {
                                            PrimaryPoles[0][k]=XX;
                                            SortFlag=FALSE;
                                    } //endif
                            } //endfor

                            if (SortFlag)
                                    m_ArrPrimaryPoles.Add(PrimaryPoles[0]);
                    } //endfor
            } //endfor
    } //endif


    for (i=0;i<RBData.seqNum;i++)
            delete	[]pAminoAcidIndeces[i];

    delete	[]pAminoAcidIndeces;


    CalculatePolesScore(&m_ArrPrimaryPoles);
    SortPolesByScore(&m_ArrPrimaryPoles);

    CorrFlag=TRUE;
    while (CorrFlag) {
            CorrFlag=FALSE;
            for (i=0;i<m_ArrPrimaryPoles.GetSize()-1;i++) {
                    j=1;
                    while ((i+j)< m_ArrPrimaryPoles.GetSize() && m_ArrPrimaryPoles[i].m_Score==m_ArrPrimaryPoles[i+j].m_Score) {
                            if (m_ArrPrimaryPoles[i]==m_ArrPrimaryPoles[i+j]) {
                                    m_ArrPrimaryPoles.RemoveAt(i);
                                    j--;
                                    CorrFlag=FALSE;
                            } //endif
                            j++;
                    } //endwhile
            } //endfor
    } //endwhile


    m_ArrStickyPoles.Copy(m_ArrPrimaryPoles);

    m_ArrPrimaryPoles.Copy(m_ArrStickyPoles);
    SortPrimaryPolesByLocation(&m_ArrPrimaryPoles);

    switch (m_NP) {
            case 0:
                    m_NumberOfPoles=RandomInt(3,7);
                    break;
            case 1:
                    m_NumberOfPoles=RandomInt(8,12);
                    break;
            case 2:
                    m_NumberOfPoles=RandomInt(13,17);
                    break;
            case 3:
                    m_NumberOfPoles=RandomInt(18,25);
                    break;
            case 4:
                    m_NumberOfPoles=RandomInt(3,7);
                    break;
            case 5:
                    m_NumberOfPoles=RandomInt(8,12);
                    break;
            case 6:
                    m_NumberOfPoles=RandomInt(13,17);
                    break;
            case 7:
                    m_NumberOfPoles=RandomInt(3,7);
                    break;
            case 8:
                    m_NumberOfPoles=RandomInt(8,12);
                    break;
    } //endswitch

    switch (m_PPSelection) {
            case 0: //Homogenous
                    size=int(m_ArrPrimaryPoles.GetSize());

                    for (i=1;i<size-1;i++) {
                            for (k=0;k<RBData.seqNum;k++) {
                                    m_ArrPrimaryPoles[i][k]=i*m_InputSequences[k].m_Seq.GetLength()/(size-1);
                            } //endfor
                    } //endfor
                    break;

            case 1: //Random Walk
                    CorrFlag=FALSE;
                    while (!CorrFlag) {
                            CorrFlag=TRUE;
                            HandleWindowsMessages();

                            for (i=0;i<m_ArrPrimaryPoles.GetSize()-1;i++) {
                                    if (m_ArrPrimaryPoles[i]>m_ArrPrimaryPoles[i+1]) {
                                            if (RandomInt(0,1)==1) {
                                                    m_ArrPrimaryPoles.RemoveAt(i);
                                            } else {
                                                    m_ArrPrimaryPoles.RemoveAt(i+1);
                                            } //endif
                                            CorrFlag=FALSE;
                                            break;
                                    } //endif
                            } //endfor
                    } //endwhile
                    break;

            case 2: //	Nearest Neighborhood"
                    CorrFlag=FALSE;
                    while (!CorrFlag) {
                            CorrFlag=TRUE;
                            HandleWindowsMessages();

                            for (i=0;i<m_ArrPrimaryPoles.GetSize()-1;i++) {
                                    if (m_ArrPrimaryPoles[i]>m_ArrPrimaryPoles[i+1]) {
                                            m_ArrPrimaryPoles.RemoveAt(i+1);
                                            CorrFlag=FALSE;
                                            i--;
                                            break;
                                    } //endif
                            } //endfor
                    } //endwhile
                    break;

            case 3: //K-Means Clustering
                    CKMeanClusterer		KmClusterer;

                    KmClusterer.m_ArrInput.Copy(m_ArrPrimaryPoles);

                    int		NPR=0;

                    for (i=0;i<m_InputSequences.GetSize();i++)
                            NPR+=m_InputSequences[i].m_Seq.GetLength();

//			NPR=int(double(NPR)/int(m_Seq.GetSize())/RandomInt(5,10)+0.5);
//			m_pRBSD->m_ExtraInfoList.AddString(_T("Km SD   :")+IntToString(NPR));
                    NPR=RandomInt(75,125)*int(m_InputSequences.GetSize())/100;
                    if (m_pRBSD!=NULL) {
                            m_pRBSD->m_ExtraInfoList.AddString(_T("Km SD F :")+IntToString(NPR));
                            m_pRBSD->m_ExtraInfoList.AddString(_T("Dim :")+IntToString(int(m_InputSequences.GetSize())));
                    } //endif

//			AfxMessageBox(IntToString(NPR));
                    i=KmClusterer.Cluster(NPR);
//			KmClusterer.ShowClusters();
                    if (m_pRBSD!=NULL) {		
                            m_pRBSD->m_ExtraInfoList.AddString(_T("Km ArrPP b4 :")+IntToString(int(m_ArrPrimaryPoles.GetSize())));
                            m_pRBSD->m_ExtraInfoList.AddString(_T("Km Iter :")+IntToString(i));
                    } //endif

                    m_ArrPrimaryPoles.RemoveAll();
                    for (i=0;i<KmClusterer.m_ArrClusters.GetSize();i++) {
//				if (KmClusterer.m_ArrClusters[i].m_ClusterRadius>1)
                                    m_ArrPrimaryPoles.Append(KmClusterer.m_ArrClusters[i].m_ArrPoles);
                    } //endfor

                    SortPrimaryPolesByLocation(&m_ArrPrimaryPoles);

                    if (m_pRBSD!=NULL)
                            m_pRBSD->m_ExtraInfoList.AddString(_T("Km ArrPP af :")+IntToString(int(m_ArrPrimaryPoles.GetSize())));

                    m_ArrStickyPoles.Copy(m_ArrPrimaryPoles);

                    size=int(m_ArrPrimaryPoles.GetSize());

                    for (i=1;i<size-1;i++) {
                            for (k=0;k<RBData.seqNum;k++) {
                                    m_ArrPrimaryPoles[i][k]=i*m_InputSequences[k].m_Seq.GetLength()/(size-1);
                            } //endfor
                    } //endfor
                    break;

    } //endswitch
    if (m_InputSequences[0].m_SS.GetLength()!=0) {
            int		HH,EE,CC;
            for (i=0;i<m_ArrStickyPoles.GetSize();i++) {
                    HH=EE=CC=XX=0;
                    for (k=0;k<RBData.seqNum;k++) {
                            j=m_ArrStickyPoles[i][k];
                            str=m_InputSequences[k].m_SS[j];
                            if (str[0]=='H')	HH++;
                            if (str[0]=='E')	EE++;
                            if (str[0]=='C')	CC++;
                            if (str[0]=='X')	XX++;
                    } //endfor

                    if (double(HH)/(RBData.seqNum-XX+1)>=0.5 || 
                            double(EE)/(RBData.seqNum-XX+1)>=0.5 ||
                            double(CC)/(RBData.seqNum-XX+1)>=0.5 ) {
                            continue;
                    } else {
                            m_ArrStickyPoles.RemoveAt(i);
                            i--;
                    } //endif
            } //endfor
    } //endif 

    while (m_ArrPrimaryPoles.GetSize()<m_NumberOfPoles) {
            XX=RandomInt(0,int(m_ArrPrimaryPoles.GetSize())-1);
            APrimaryPole=m_ArrPrimaryPoles.GetAt(XX);
            m_ArrPrimaryPoles.InsertAt(XX,APrimaryPole);
    } //endwhile

    while (m_ArrPrimaryPoles.GetSize()>m_NumberOfPoles) {
            XX=RandomInt(0,int(m_ArrPrimaryPoles.GetSize())-1);
            m_ArrPrimaryPoles.RemoveAt(XX);
    } //endwhile

    for (i=0;i<RBData.seqNum;i++)
            APrimaryPole[i]=0;
    m_ArrPrimaryPoles.InsertAt(0,APrimaryPole);
    m_ArrStickyPoles.InsertAt(0,APrimaryPole);

    for (i=0;i<RBData.seqNum;i++)
            APrimaryPole[i]=m_InputSequences[i].m_Seq.GetLength();

    m_ArrPrimaryPoles.Add(APrimaryPole);
    m_ArrStickyPoles.Add(APrimaryPole);

    if (pCounters != NULL)
        free (pCounters);
    pCounters = NULL;
    if (pCounterLimits != NULL)
        free (pCounterLimits);
    pCounterLimits = NULL;
    if (PrimaryPoles != NULL) {
        free (PrimaryPoles[0]);
        free (PrimaryPoles[1]);
        free (PrimaryPoles[2]);
        free (PrimaryPoles[3]);
        free (PrimaryPoles[4]);
        free (PrimaryPoles);
    }
    PrimaryPoles = NULL;
    if (APrimaryPole != NULL)
        free (APrimaryPole);
    APrimaryPole = NULL;
    
}

void FixPrimaryPoles(RBAlignData RBData) {
    MOATypeDimn k;
    MOATypeShape i;
    int CorrFlag;

    CorrFlag= 0;
    while (CorrFlag == 0) {
        CorrFlag=1;
        for (i=1;i<RBData.pArrPrimaryPolesNum-1;i++) {
            for (k=0;k<pData->seqNum;k++) {
                if (RBData.pArrPrimaryPoles[i][k]>RBData.pArrPrimaryPoles[i+1][k]) {
                    switch (RBData.FixType) {
                        case 0:
                            RBData.pArrPrimaryPoles[i][k]=(RBData.pArrPrimaryPoles[i-1][k]+RBData.pArrPrimaryPoles[i+1][k])/2;
                            break;
                        case 1:
                            RBData.pArrPrimaryPoles[i][k]=RBData.pArrPrimaryPoles[i+1][k];
                            break;
                    } //endswtich
                    CorrFlag=0;
                    break;
                } //endif
            } //endfor
        } //endfor
    } //endwhile
}

void RBSolve(CArray <SProteinSequence,SProteinSequence&> *pSeqs, CString ScoreMatrix, int GapOpening, int GapExtension, int ScoreMode) {
    int	i,j,score,NumOfTries;
    BOOL flag;
    CString str;
    CTime StartTime,FinishTime,StartAllTime,FinishAllTime;
    CTimeSpan ElapsedTime;
    CPole DummyPole;
    int	NextMNP,NumberOfRounds;
    BOOL ChangeFlag;
    CArray <SProteinSequence,SProteinSequence&>	BestAnsSeq;
    double BestBALiBASEScore;
    CArray <int,int> ArrScore;
    CArray <double,double> ArrPCRes,ArrTCScore,ArrMSScore;

    StartAllTime=CTime::GetCurrentTime();

    m_PPGeneration=2;
    m_PPSelection=0;
    m_NP=2;
    m_PPWideness=1;
    m_PPSimilarity=3;
    m_RunningMode=4;
    m_TryMode	=2;

    m_Dimension=int(pSeqs->GetSize());
    m_SticknessPercentage=100;

    switch (m_RunningMode) {
        case 0:
        case 1:
            m_NumOfRuns=1;
            break;
        case 2:
            m_NumOfRuns=5;
            break;
        case 3:
            m_NumOfRuns=10;
        case 4:
            m_NumOfRuns=10;
    } //endswitch

    m_InputSequences.Copy(*pSeqs);
    m_AlignedSequences.SetSize(m_InputSequences.GetSize());
    for (i=0;i<m_InputSequences.GetSize();i++)
        m_AlignedSequences[i].m_Name=m_InputSequences[i].m_Name;

    switch (m_NP) {
        case 0:
        case 1:
        case 2:
            NumberOfRounds=2;
            break;
        case 3:
            NumberOfRounds=1;
            break;
        case 4:
        case 5:
        case 6:
            NumberOfRounds=2;
            break;
        case 7:
        case 8:
            NumberOfRounds=3;
            break;
    } //endswitch

    m_ASC.Initialize(ScoreMatrix,GapOpening,GapExtension);
    m_ASC.CalculateWeightMatrix(ScoreMode,&m_InputSequences); /*do DP pairwise alignment for each pair of sequences and calculate the pair-wise score weight*/

    m_TryCounter=0;
    for (i=0;i<m_NumOfRuns;i++) {
        StartTime=CTime::GetCurrentTime();
        FindInitialPrimaryPoles();
        switch (m_TryMode) {
            case 0:
                m_NumberOfTries=5;
                break;
            case 1:
                m_NumberOfTries=10;
                break;
            case 2:
                m_NumberOfTries=15;
                break;
            case 3:
                m_NumberOfTries=int(m_NumberOfPoles/2)+1;
                break;
            case 4:
                m_NumberOfTries=m_NumberOfPoles;
                break;
        } //endswitch
        FixPrimaryPoles(&m_ArrPrimaryPoles);

        m_ArrPoles.Copy(m_ArrPrimaryPoles);
        NumOfTries=int(1+double(3*m_TryCounter)/m_NumberOfTries);
        FindConnectingPoles(&m_ArrPoles,NumOfTries);
        MakeSequencesFromPoles(&m_AlignedSequences,&m_ArrPoles);
        SeperatePrimaryPoles();

            m_RoundCounter=0;
            while (m_RoundCounter<NumberOfRounds) {
                    m_TryCounter=0;
                    ChangeFlag=TRUE;

                    while (m_TryCounter<(m_NumberOfTries*(m_RoundCounter+1)/NumberOfRounds)) {
                            ChangeFlag=FALSE;
                            if (MoveBlocks())				ChangeFlag=TRUE;
                            if (MovePrimaryPoles())			ChangeFlag=TRUE;
                            if (JumpToStickyPoles())		ChangeFlag=TRUE;
                            if (StickPrimaryPolesToEachOther())			ChangeFlag=TRUE;
                            if (AlignPrimaryPoles())		ChangeFlag=TRUE;

                            if (ChangeFlag)	m_TryCounter=0;
                            else			m_TryCounter++;

                            if (m_pRBSD!=NULL)
                                    m_pRBSD->m_TryCounterProgress.SetPos(m_TryCounter);

                            if (m_TryCounter!=0) {
                                    str =_T("    R ");
                                    str+=IntToString(m_RoundCounter+1);
                                    str+=_T(" - Try : ");
                                    str+=IntToString(m_TryCounter);
                                    str+=_T(" / ");
                                    str+=IntToString(m_NumberOfTries-1);
                                    str+=_T(" : ");
                                    str+=IntToString(int(1+double(3*m_TryCounter)/m_NumberOfTries));
                                    if (m_pRBSD!=NULL) {
                                            m_pRBSD->m_EnergyList.AddString(str);
                                            m_pRBSD->m_EnergyList.SetCurSel(m_pRBSD->m_EnergyList.GetCount()-1);
                                    } //endif
                            } //endif

                            FinishTime=CTime::GetCurrentTime();
                            ElapsedTime=FinishTime-StartTime;
                            str=ElapsedTime.Format("%H h, %M m, %S s" );
                            if (m_pRBSD!=NULL) {
                                    m_pRBSD->m_csElapsedTime=str;
                                    m_pRBSD->UpdateData(FALSE);
                            } //endif

                            if (m_RunningMode==0)
                                    if (IDNO==AfxMessageBox(_T("Do you want to continue "),MB_YESNO|MB_ICONQUESTION)) {
                                            AssignSecondaryStructure();
                                            return;
                                    } //endif

                    } //endwhile

                    switch (m_NP) {
                            case 0:
                            case 1:
                            case 2:
                            case 3:
                                    NextMNP=0;
                                    break;
                            case 4:
                                    if (m_RoundCounter==0)	NextMNP=RandomInt(8,12);
                                    break;
                            case 5:
                                    if (m_RoundCounter==0)	NextMNP=RandomInt(13,17);
                                    break;
                            case 6:
                                    if (m_RoundCounter==0)	NextMNP=RandomInt(18,25);
                                    break;
                            case 7:
                                    if (m_RoundCounter==0)	NextMNP=RandomInt(8,12);
                                    if (m_RoundCounter==1)	NextMNP=RandomInt(13,17);
                                    break;
                            case 8:
                                    if (m_RoundCounter==1)	NextMNP=RandomInt(13,17);
                                    if (m_RoundCounter==0)	NextMNP=RandomInt(18,25);
                                    break;

                    } //endswitch

                    while (m_ArrPrimaryPoles.GetSize()<NextMNP) {
                            j=RandomInt(1,int(m_ArrPrimaryPoles.GetSize())-2);
                            DummyPole=m_ArrPrimaryPoles[j];
                            m_ArrPrimaryPoles.InsertAt(j,DummyPole);
                    } //endwhile

                    m_RoundCounter++;
                    if (m_RoundCounter<NumberOfRounds) {
                            m_NumberOfPoles=NextMNP;
                            if (m_pRBSD!=NULL) {
                                    m_pRBSD->m_ExtraInfoList.AddString(_T("NP :")+IntToString(int(m_ArrPrimaryPoles.GetSize())));
                                    m_pRBSD->m_ExtraInfoList.SetCurSel(m_pRBSD->m_ExtraInfoList.GetCount()-1);
                            } //endif
                    } //endif
            } //endwhile

            CalculatePolesScore(&m_ArrPrimaryPoles);
//		ShowPolesScore(&m_ArrPrimaryPoles);

//		if (EliminateUnnecassaryPrimaryPoles())		ChangeFlag=TRUE;
            if (EliminateUnnecassaryPoles())	ChangeFlag=TRUE;

            MakeSequencesFromPoles(&m_AlignedSequences,&m_ArrPoles);
            CalculateBALiBASEScore(&m_AlignedSequences);

            flag=FALSE;
            if (i==0)						flag=TRUE;
            else if (m_BBS.m_PCRes>BestBALiBASEScore)	flag=TRUE;

            if (flag) {
                    BestBALiBASEScore=m_BBS.m_PCRes;
                    BestAnsSeq.Copy(m_AlignedSequences);
            } //endif

            score=CalculateAlignmentScore(&m_AlignedSequences);
            ArrScore.Add(score);
            ArrPCRes.Add(m_BBS.m_PCRes);
            ArrTCScore.Add(m_BBS.m_TCScore);
            ArrMSScore.Add(m_BBS.m_MSScore);

            if (m_pRBSD!=NULL) {		
                    m_pRBSD->m_ExtraInfoList.AddString(IntToString(i)+_T(" -->  ")+IntToString(score));
                    m_pRBSD->m_ExtraInfoList.SetCurSel(m_pRBSD->m_ExtraInfoList.GetCount()-1);
            } //endif


            flag=FALSE;
            while (!flag) {
                    flag=TRUE;

                    for (j=0;j<m_ArrPrimaryPoles.GetSize()-1;j++)
                            if (m_ArrPrimaryPoles[j]==m_ArrPrimaryPoles[j+1]) {
                                    m_ArrPrimaryPoles.RemoveAt(j+1);
                                    flag=FALSE;
                            } //endif
            } //endwhile

            if (m_pRBSD!=NULL) {
                    m_pRBSD->m_ExtraInfoList.AddString(_T("PP After Optimization"));
                    m_pRBSD->m_ExtraInfoList.AddString(_T("NP :")+IntToString(int(m_ArrPrimaryPoles.GetSize())));
            } //endif

            FinishTime=CTime::GetCurrentTime();
            ElapsedTime=FinishTime-StartTime;
            str=ElapsedTime.Format("%H h, %M m, %S s" );

            if (m_pRBSD!=NULL) {
                    m_pRBSD->m_ExtraInfoList.AddString(_T("Optimization Time"));
                    m_pRBSD->m_ExtraInfoList.AddString(str);
                    m_pRBSD->m_ExtraInfoList.AddString(_T("       "));
                    m_pRBSD->m_EnergyList.AddString(_T("       "));
                    m_pRBSD->m_ExtraInfoList.SetCurSel(m_pRBSD->m_ExtraInfoList.GetCount()-1);
            } //endif
    } //endfor


    FinishAllTime=CTime::GetCurrentTime();
    ElapsedTime=FinishAllTime-StartAllTime;
    str=ElapsedTime.Format("%H h, %M m, %S s" );

    if (m_pRBSD!=NULL) {
            m_pRBSD->m_ExtraInfoList.AddString(_T("       "));
            m_pRBSD->m_ExtraInfoList.AddString(_T("       "));
            m_pRBSD->m_ExtraInfoList.AddString(_T("Optimization All Time"));
            m_pRBSD->m_ExtraInfoList.AddString(str);
            m_pRBSD->m_ExtraInfoList.SetCurSel(m_pRBSD->m_ExtraInfoList.GetCount()-1);
    } //endif


    double	dd1,dd2,dd3,dd4,MaxSPS,MinSPS;
    int		MaxScore;

    dd1=dd2=dd3=dd4=0;
    MaxScore=ArrScore[0];
    MaxSPS=MinSPS=ArrPCRes[0];

    for (i=1;i<ArrScore.GetSize();i++) {
            if (MinSPS>ArrPCRes[i])			MinSPS=ArrPCRes[i];
            if (MaxSPS<ArrPCRes[i])			MaxSPS=ArrPCRes[i];
            if (MaxScore<ArrScore[i])		MaxScore=ArrScore[i];
    } //endfor

    str=_T("");
    for (i=0;i<ArrScore.GetSize();i++) {
            str+=IntToString(ArrScore[i]);
            if (int(ArrScore[i])==MaxScore)	str+='*';
            str+=_T("		");

            str+=DoubleToString(ArrPCRes[i]);
            if (ArrPCRes[i]==MaxSPS)	str+=_T(" ++");
            if (ArrPCRes[i]==MinSPS)	str+=_T(" ---");
//		str+=_T("	");

//		str+=DoubleToString(ArrTCScore[i]);
//		str+=_T("	");
//		str+=DoubleToString(ArrMSScore[i]);
            str+=_T("\n");

            dd1+=ArrScore[i];
            dd2+=ArrPCRes[i];
            dd3+=ArrTCScore[i];
            dd4+=ArrMSScore[i];
    } //endfor

    dd1/=ArrScore.GetSize();
    dd2/=ArrScore.GetSize();
    dd3/=ArrScore.GetSize();
    dd4/=ArrScore.GetSize();
    str+=_T("\n");
    str+=_T("----------------------------------\n");
    str+=IntToString(int(dd1));
    str+=_T("		");
    str+=DoubleToString(dd2);
//	str+=_T("	");
//	str+=DoubleToString(dd3);
//	str+=_T("	");
//	str+=DoubleToString(dd4);



    m_AlignedSequences.Copy(BestAnsSeq);
    AssignSecondaryStructure();
}