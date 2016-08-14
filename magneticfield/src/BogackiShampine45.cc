#include "BogackiShampine45.hh"
#include "G4LineSection.hh"
#include "G4ClassicalRK4.hh"

#include "G4SystemOfUnits.hh"

#define ncomp G4FieldTrack::ncompSVEC

//Constructor
BogackiShampine45::BogackiShampine45(G4EquationOfMotion *EqRhs,
                                     G4int     noIntegrationVariables)
   : G4MagIntegratorStepper(EqRhs, std::min(noIntegrationVariables,G4int(ncomp))),
     fInterpolationPrepared(false)
{
    Init();
}

void BogackiShampine45::Init()
{
/*
    This pair is from "An Efficient Runge-Kutta (4,5) Pair" by P. Bogacki
    and L.F. Shampine, Rept. 89-20, Math. Dept., Southern Methodist
    University, Dallas, Texas, USA, 1989.  The authors are grateful to
    P. Bogacki for his assistance in implementing the pair.  Shampine and
    Bogacki subsequently modified the formula to enhance the reliability of
    the pair.  The original fourth order formula is used in an estimate of
    the local error.  If the step fails, the computation is broken off.  If
    the step is acceptable, the first evaluation of the next step is done,
    i.e., the pair is implemented as FSAL and the local error of the step
    is again estimated with a fourth order formula using the additional data.
    The step must succeed with both estimators to be accepted.  When the
    second estimate is formed, it is used for the subsequent adjustment of
    the step size because it is of higher quality.  The two fourth order
    formulas are well matched to leading order, and only exceptionally do
    the estimators disagree -- problems with discontinuous coefficients are
    handled more reliably by using two estimators as is global error
    estimation.
*/
    /*
          nstage = 8;
          rkcom5.fsal = true;
          rkcom5.order = 4;
          rkcom5.tanang = 5.2e0;
          rkcom5.stbrad = 3.9e0;
          rkcom5.safety = 0.8e0;
          intp = true;
          reqstg = true;
          mintp = 6;
          lintpl = 6;
          rkcom5.nsec = 2;

          ptr[1] = 0;
          ptr[2] = 1;
          ptr[3] = 2;
          ptr[4] = 3;
          ptr[5] = 4;
          ptr[6] = 5;
          ptr[7] = 6;
          ptr[8] = 7;
*/
          a[2][1] = 1.0e0/6.0e0;
          a[3][1] = 2.e0/27.e0;
          a[3][2] = 4.e0/27.e0;
          a[4][1] = 183.e0/1372.e0;
          a[4][2] = -162.e0/343.e0;
          a[4][3] = 1053.e0/1372.e0;
          a[5][1] = 68.e0/297.e0;
          a[5][2] = -4.e0/11.e0;
          a[5][3] = 42.e0/143.e0;
          a[5][4] = 1960.e0/3861.e0;
          a[6][1] = 597.e0/22528.e0;
          a[6][2] = 81.e0/352.e0;
          a[6][3] = 63099.e0/585728.e0;
          a[6][4] = 58653.e0/366080.e0;
          a[6][5] = 4617.e0/20480.e0;
          a[7][1] = 174197.e0/959244.e0;
          a[7][2] = -30942.e0/79937.e0;
          a[7][3] = 8152137.e0/19744439.e0;
          a[7][4] = 666106.e0/1039181.e0;
          a[7][5] = -29421.e0/29068.e0;
          a[7][6] = 482048.e0/414219.e0;
          a[8][1] = 587.e0/8064.e0;
          a[8][2] = 0.e0;
          a[8][3] = 4440339.e0/15491840.e0;
          a[8][4] = 24353.e0/124800.e0;
          a[8][5] = 387.e0/44800.e0;
          a[8][6] = 2152.e0/5985.e0;
          a[8][7] = 7267.e0/94080.e0;
/*
     The coefficients B(*) refer to the formula of order 4.
*/
          b[1] = 2479.e0/34992.e0;
          b[2] = 0.e0;
          b[3] = 123.e0/416.e0;
          b[4] = 612941.e0/3411720.e0;
          b[5] = 43.e0/1440.e0;
          b[6] = 2272.e0/6561.e0;
          b[7] = 79937.e0/1113912.e0;
          b[8] = 3293.e0/556956.e0;
/*
      The coefficients E(*) refer to an estimate of the local error based on
      the first formula of order 4.  It is the difference of the fifth order
      result, here located in A(8,*), and the fourth order result.  By
      construction both E(2) and E(7) are zero.
*/
          e[1] = -3.e0/1280.e0;
          e[2] = 0.e0;
          e[3] = 6561.e0/632320.e0;
          e[4] = -343.e0/20800.e0;
          e[5] = 243.e0/12800.e0;
          e[6] = -1.e0/95.e0;
          e[7] = 0.e0;
/*
          c[1] = 0.e0;
          c[2] = 1.e0/6.e0;
          c[3] = 2.e0/9.e0;
          c[4] = 3.e0/7.e0;
          c[5] = 2.e0/3.e0;
          c[6] = 3.e0/4.e0;
          c[7] = 1.e0;
          c[8] = 1.e0;
*/
/*
      To do interpolation with this pair, some extra stages have to be computed.
      The following additional A(*,*) and C(*) coefficients are for this purpose.
      In addition there is an array R(*,*) that plays a role for interpolation
      analogous to that of BHAT(*) for the basic step.
*/
/*
          c[9] = 1.e0/2.e0;
          c[10] = 5.e0/6.e0;
          c[11] = 1.e0/9.e0;
*/
          a[9][1] = 455.e0/6144.e0;
          a[10][1] = -837888343715.e0/13176988637184.e0;
          a[11][1] = 98719073263.e0/1551965184000.e0;
          a[9][2] = 0.e0;
          a[10][2] = 30409415.e0/52955362.e0;
          a[11][2] = 1307.e0/123552.e0;
          a[9][3] = 10256301.e0/35409920.e0;
          a[10][3] = -48321525963.e0/759168069632.e0;
          a[11][3] = 4632066559387.e0/70181753241600.e0;
          a[9][4] = 2307361.e0/17971200.e0;
          a[10][4] = 8530738453321.e0/197654829557760.e0;
          a[11][4] = 7828594302389.e0/382182512025600.e0;
          a[9][5] = -387.e0/102400.e0;
          a[10][5] = 1361640523001.e0/1626788720640.e0;
          a[11][5] = 40763687.e0/11070259200.e0;
          a[9][6] = 73.e0/5130.e0;
          a[10][6] = -13143060689.e0/38604458898.e0;
          a[11][6] = 34872732407.e0/224610586200.e0;
          a[9][7] = -7267.e0/215040.e0;
          a[10][7] = 18700221969.e0/379584034816.e0;
          a[11][7] = -2561897.e0/30105600.e0;
          a[9][8] = 1.e0/32.e0;
          a[10][8] = -5831595.e0/847285792.e0;
          a[11][8] = 1.e0/10.e0;
          a[10][9] = -5183640.e0/26477681.e0;
          a[11][9] = -1.e0/10.e0;
          a[11][10] = -1403317093.e0/11371610250.e0;

          for (G4int i = 1; i <= 11; i++) {
             r[i][1] = 0.e0;
          }
          for (G4int i = 1; i <= 6; i++) {
             r[2][i] = 0.e0;
          }
          r[1][6] = -12134338393.e0/1050809760.e0;
          r[1][5] = -1620741229.e0/50038560.e0;
          r[1][4] = -2048058893.e0/59875200.e0;
          r[1][3] = -87098480009.e0/5254048800.e0;
          r[1][2] = -11513270273.e0/3502699200.e0;

          r[3][6] = -33197340367.e0/1218433216.e0;
          r[3][5] = -539868024987.e0/6092166080.e0;
          r[3][4] = -39991188681.e0/374902528.e0;
          r[3][3] = -69509738227.e0/1218433216.e0;
          r[3][2] = -29327744613.e0/2436866432.e0;

          r[4][6] = -284800997201.e0/19905339168.e0;
          r[4][5] = -7896875450471.e0/165877826400.e0;
          r[4][4] = -333945812879.e0/5671036800.e0;
          r[4][3] = -16209923456237.e0/497633479200.e0;
          r[4][2] = -2382590741699.e0/331755652800.e0;

          r[5][6] = -540919.e0/741312.e0;
          r[5][5] = -103626067.e0/43243200.e0;
          r[5][4] = -633779.e0/211200.e0;
          r[5][3] = -32406787.e0/18532800.e0;
          r[5][2] = -36591193.e0/86486400.e0;

          r[6][6] = 7157998304.e0/374350977.e0;
          r[6][5] = 30405842464.e0/623918295.e0;
          r[6][4] = 183022264.e0/5332635.e0;
          r[6][3] = -3357024032.e0/1871754885.e0;
          r[6][2] = -611586736.e0/89131185.e0;

          r[7][6] = -138073.e0/9408.e0;
          r[7][5] = -719433.e0/15680.e0;
          r[7][4] = -1620541.e0/31360.e0;
          r[7][3] = -385151.e0/15680.e0;
          r[7][2] = -65403.e0/15680.e0;

          r[8][6] = 1245.e0/64.e0;
          r[8][5] = 3991.e0/64.e0;
          r[8][4] = 4715.e0/64.e0;
          r[8][3] = 2501.e0/64.e0;
          r[8][2] = 149.e0/16.e0;
          r[8][1] = 1.e0;

          r[9][6] = 55.e0/3.e0;
          r[9][5] = 71.e0;
          r[9][4] = 103.e0;
          r[9][3] = 199.e0/3.e0;
          r[9][2] = 16.0e0;

          r[10][6] = -1774004627.e0/75810735.e0;
          r[10][5] = -1774004627.e0/25270245.e0;
          r[10][4] = -26477681.e0/359975.e0;
          r[10][3] = -11411880511.e0/379053675.e0;
          r[10][2] = -423642896.e0/126351225.e0;

          r[11][6] = 35.e0;
          r[11][5] = 105.e0;
          r[11][4] = 117.e0;
          r[11][3] = 59.e0;
          r[11][2] = 12.e0;
}

BogackiShampine45::~BogackiShampine45()
{
}


void BogackiShampine45::Stepper(const G4double yInput[],
                                const G4double dydx[],
                                      G4double Step,
                                      G4double yOutput[],
                                      G4double yErr[])
{

    memcpy(yIn,yInput,sizeof(G4double)*ncomp);
    memcpy(dydxIn,dydx,sizeof(G4double)*ncomp);

    //copy time to temp
    yTemp[7] = yIn[7];


    // RightHandSide(yIn, dydx) ;
    // 1st Step - Not doing, getting passed

    for(G4int i=0; i<GetNumberOfVariables(); ++i)
    {
        yTemp[i] = yIn[i] + a[2][1]*Step*dydx[i] ;
    }
    RightHandSide(yTemp, ak2);              // 2nd Step

    for(G4int i=0; i<GetNumberOfVariables(); ++i)
    {
        yTemp[i] = yIn[i] + Step*(a[3][1]*dydx[i] + a[3][2]*ak2[i]) ;
    }
    RightHandSide(yTemp, ak3);              // 3rd Step

    for(G4int i=0; i<GetNumberOfVariables(); ++i)
    {
        yTemp[i] = yIn[i] + Step*(a[4][1]*dydx[i] + a[4][2]*ak2[i] + a[4][3]*ak3[i]) ;
    }
    RightHandSide(yTemp, ak4);              // 4th Step

    for(G4int i=0; i<GetNumberOfVariables(); ++i)
    {
        yTemp[i] = yIn[i] + Step*(a[5][1]*dydx[i] + a[5][2]*ak2[i] + a[5][3]*ak3[i] +
                                  a[5][4]*ak4[i]) ;
    }
    RightHandSide(yTemp, ak5);              // 5th Step

    for(G4int i=0; i<GetNumberOfVariables(); ++i)
    {
        yTemp[i] = yIn[i] + Step*(a[6][1]*dydx[i] + a[6][2]*ak2[i] + a[6][3]*ak3[i] +
                                  a[6][4]*ak4[i] + a[6][5]*ak5[i]) ;
    }
    RightHandSide(yTemp, ak6);              // 6th Step

    for(G4int i=0; i<GetNumberOfVariables(); ++i)
    {
        yTemp[i] = yIn[i] + Step*(a[7][1]*dydx[i] + a[7][2]*ak2[i] + a[7][3]*ak3[i] +
                                  a[7][4]*ak4[i] + a[7][5]*ak5[i] + a[7][6]*ak6[i]);
    }
    RightHandSide(yTemp, ak7);				//7th Step

    for(G4int i=0; i<GetNumberOfVariables(); ++i)
    {
        yOutput[i] = yIn[i] + Step*(a[8][1]*dydx[i] + a[8][2]*ak2[i] + a[8][3]*ak3[i] +
                                 a[8][4]*ak4[i] + a[8][5]*ak5[i] + a[8][6]*ak6[i] +
                                 a[8][7]*ak7[i]);
    }
    RightHandSide(yOutput, ak8);				//8th Step - Final one Using FSAL

    memcpy(dydxOut,ak8,sizeof(G4double)*G4FieldTrack::ncompSVEC);
    memcpy(yOut,yOutput,sizeof(G4double)*G4FieldTrack::ncompSVEC);

    for(G4int i=0; i<GetNumberOfVariables(); ++i)
    {
        yErr[i] = Step*(e[1]*dydx[i] + e[3]*ak3[i] + e[4]*ak4[i] +
                        e[5]*ak5[i] + e[6]*ak6[i]);
    }

    memcpy(yError,yErr,sizeof(G4double)*ncomp);

    fLastStepLength = Step;
    fInterpolationPrepared = false;
}

G4double  BogackiShampine45::DistChord() const
{

    // Using interpolation, requires only 3 extra stages (ie 3 extra field evaluations )
    Interpolate(fLastStepLength/2., yMid);
/*
    //check the interpolation!
    G4double tmp[12],err[12];
    BogackiShampine45* This = const_cast<BogackiShampine45*>(this);
    G4ClassicalRK4 classicalRK(This->GetEquationOfMotion());
    classicalRK.Stepper(yIn,dydxIn,fLastStepLength/2.,tmp,err);

    for (G4int i = 0; i < GetNumberOfVariables(); ++i)
    {
        tmp[i] -= yMid[i];
    }

    G4double error = sqr(tmp[0]) + sqr(tmp[1]) + sqr(tmp[2]);
    error /= sqr(fLastStepLength);

    G4double estError = sqr(yError[0]) + sqr(yError[1]) + sqr(yError[2]);
    estError /= sqr(fLastStepLength);

    if (error > 1e-5)
        G4cout<<"error: "<<error<<" estError: "<<estError<<G4endl;
*/

    // Use stored values of Initial and Endpoint + new Midpoint to evaluate  distance of Chord
    G4ThreeVector initialPoint(yIn[0], yIn[1], yIn[2]);
    G4ThreeVector finalPoint(yOut[0], yOut[1], yOut[2]);
    G4ThreeVector midPoint(yMid[0], yMid[1], yMid[2]);

    G4double distChord;
    if (initialPoint != finalPoint)
    {
        distChord = G4LineSection::Distline( midPoint, initialPoint, finalPoint );
    }
    else
    {
        distChord = (midPoint-initialPoint).mag();
    }
    //G4cout<<"distChord: "<<distChord<<G4endl;
    return distChord;
}

void BogackiShampine45::FormInterpolation()
{
    yTemp[7]  = yIn[7];

    //Evaluate the extra stages :
    for(G4int i=0; i<GetNumberOfVariables(); ++i)
    {
        yTemp[i] = yIn[i] + fLastStepLength*(a[9][1]*dydxIn[i] + a[9][2]*ak2[i] + a[9][3]*ak3[i] +
                                  a[9][4]*ak4[i] + a[9][5]*ak5[i] + a[9][6]*ak6[i] +
                                  a[9][7]*ak7[i] + a[9][8]*ak8[i]);
    }

    RightHandSide(yTemp, ak9);		//9th stage

    for(G4int i=0; i<GetNumberOfVariables(); ++i)
    {
        yTemp[i] = yIn[i] + fLastStepLength*(a[10][1]*dydxIn[i] + a[10][2]*ak2[i] + a[10][3]*ak3[i] +
                                  a[10][4]*ak4[i] + a[10][5]*ak5[i] + a[10][6]*ak6[i] +
                                  a[10][7]*ak7[i] + a[10][8]*ak8[i] + a[10][9]*ak9[i]);
    }

    RightHandSide(yTemp, ak10);		//10th stage

    for(G4int i=0; i<GetNumberOfVariables(); ++i)
    {
        yTemp[i] = yIn[i] + fLastStepLength*(a[11][1]*dydxIn[i] + a[11][2]*ak2[i] + a[11][3]*ak3[i] +
                                  a[11][4]*ak4[i] + a[11][5]*ak5[i] + a[11][6]*ak6[i] +
                                  a[11][7]*ak7[i] + a[11][8]*ak8[i] + a[11][9]*ak9[i] +
                                  a[11][10]*ak10[i]);
    }
    RightHandSide(yTemp, ak11);		//11th stage

    //  In future we can restrict the number of variables interpolated
    G4int nwant = GetNumberOfVariables();

//  Form the coefficients of the interpolating polynomial in its shifted
//  and scaled form.  The terms are grouped to minimize the errors
//  of the transformation, to cope with ill-conditioning. ( From RKSUITE )
//
    for (G4int l = 0; l < nwant; ++l)
    {
        //  Coefficient of tau^6
        p[4][l] =   r[5][6]*ak5[l] +
                  ((r[10][6]*ak10[l] + r[8][6]*ak8[l]) +
                   (r[7][6]*ak7[l] + r[6][6]*ak6[l]))  +
                  ((r[4][6]*ak4[l] + r[9][6]*ak9[l]) +
                   (r[3][6]*ak3[l] + r[11][6]*ak11[l]) +
                    r[1][6]*dydxIn[l]);
        //  Coefficient of tau^5
        p[3][l] = (r[10][5]*ak10[l] + r[9][5]*ak9[l])  +
                 ((r[7][5]*ak7[l] + r[6][5]*ak6[l]) +
                   r[5][5]*ak5[l])  +  ((r[4][5]*ak4[l] +
                   r[8][5]*ak8[l]) + (r[3][5]*ak3[l] +
                   r[11][5]*ak11[l]) + r[1][5]*dydxIn[l]);
        //  Coefficient of tau^4
        p[2][l] = ((r[4][4]*ak4[l] + r[8][4]*ak8[l]) +
                   (r[7][4]*ak7[l] + r[6][4]*ak6[l]) +
                    r[5][4]*ak5[l]) + ((r[10][4]*ak10[l] +
                    r[9][4]*ak9[l]) +  (r[3][4]*ak3[l] +
                    r[11][4]*ak11[l]) + r[1][4]*dydxIn[l]);
        //  Coefficient of tau^3
        p[1][l] =  r[5][3]*ak5[l] + r[6][3]*ak6[l]  +
                 ((r[3][3]*ak3[l] + r[9][3]*ak9[l]) +
                 (r[10][3]*ak10[l]+ r[8][3]*ak8[l]) + r[1][3]*dydxIn[l]) +
                 ((r[4][3]*ak4[l] + r[11][3]*ak11[l]) + r[7][3]*ak7[l]);
        //  Coefficient of tau^2
        p[0][l] = r[5][2]*ak5[l]  + ((r[6][2]*ak6[l] +
                  r[8][2]*ak8[l]) +   r[1][2]*dydxIn[l])  +
                ((r[3][2]*ak3[l]  +   r[9][2]*ak9[l]) +
                 r[10][2]*ak10[l])+ ((r[4][2]*ak4[l] +
                 r[11][2]*ak11[l]) +   r[7][2]*ak7[l]);
      }

      //Scale all the coefficients by the step size.
      for (G4int i = 0; i < 5; ++i)
      {
         for (G4int icomp = 0; icomp < nwant; ++icomp)
         {
            p[i][icomp] *= fLastStepLength;
         }
      }

    fInterpolationPrepared = true;
}



void BogackiShampine45::Interpolate(G4double step, G4double yOutput[]) const
{
    if(!fInterpolationPrepared)
    {
        //call non-const function inside const function
        const_cast<BogackiShampine45*>(this)->FormInterpolation();
    }

    //we use end point as a start point for interpolation!
    G4double tau = step/fLastStepLength - 1.;

    if (tau < -1-perMillion || tau > perMillion)
    {
        char buff[256];
        sprintf(buff, "tau(%g) is out of interpolation interval(-1,0)!",tau);
        G4Exception("BogackiShampine45::Interpolate()", "GeomField0001",
              FatalException, buff);

    }
    const G4int nwant = GetNumberOfVariables();
    const G4int norder = 6;

    for (G4int icomp = 0; icomp < nwant; ++icomp)
    {
      yOutput[icomp] = p[norder-2][icomp] * tau;
    }
    for (G4int k = norder - 2; k >= 1; --k)
    {
      for (G4int icomp = 0; icomp < nwant; ++icomp)
      {
         yOutput[icomp] = (yOutput[icomp] + p[k-1][icomp]) * tau;
      }
    }
    for (G4int icomp = 0; icomp < nwant; ++icomp)
    {
      yOutput[icomp] = (yOutput[icomp] + fLastStepLength * dydxOut[icomp]) * tau + yOut[icomp];
    }
}
