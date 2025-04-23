/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::arTransport<Thermo>::arTransport
(
    const Thermo& t,
    const scalar mu,
    const scalar Pr
)
:
    Thermo(t),
    mu_(mu),
    rPr_(1.0/Pr)
{}


template<class Thermo>
inline Foam::arTransport<Thermo>::arTransport
(
    const word& name,
    const arTransport& ct
)
:
    Thermo(name, ct),
    mu_(ct.mu_),
    rPr_(ct.rPr_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::arTransport<Thermo> >
Foam::arTransport<Thermo>::clone() const
{
    return autoPtr<arTransport<Thermo> >
    (
        new arTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::arTransport<Thermo> >
Foam::arTransport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<arTransport<Thermo> >
    (
        new arTransport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::arTransport<Thermo> >
Foam::arTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<arTransport<Thermo> >
    (
        new arTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::arTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return mu_;
}


template<class Thermo>
inline Foam::scalar Foam::arTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
	//old version
    //return this->Cp(p, T)*mu(p, T)*rPr_;
    // Thermal conductivity kappa [W/(m.K)] function of T, for Argon plasma, 
   // tabulated for T from T0=200K to 20000K with tabulation interval of dT=100K

    int i_index;
    scalar dT=100;
    scalar T0=200;
    scalar Temp_Argon;
    scalar kappa_T_Argon;

   // check that T is within the tabulation range
    if (T> 20000 || T<200)
      {
          FatalErrorIn
          (
              " Attempt to use ConstTransportI,"
              " out of temperature range "
          )   << " T_kappa="<<T<<"; "
              << abort(FatalError);
      }

   // Thermal conductivity kappa [W/(m.K)] , Argon :
   // data table provided by Aubreton and Co.

    scalar kappa_Argon[199]=
       {
          .14000E-01, .18800E-01, .23200E-01, .26498E-01, .30286E-01,
          .33812E-01, .37131E-01, .40284E-01, .43298E-01, .46197E-01,
          .48995E-01, .51708E-01, .54344E-01, .56913E-01, .59421E-01,
          .61875E-01, .64279E-01, .66638E-01, .68955E-01, .71234E-01,
          .73477E-01, .75686E-01, .77864E-01, .80013E-01, .82135E-01,
          .84230E-01, .86301E-01, .88349E-01, .90375E-01, .92379E-01,
          .94363E-01, .96329E-01, .98276E-01, .10020E+00, .10212E+00,
          .10401E+00, .10589E+00, .10776E+00, .10961E+00, .11145E+00,
          .11328E+00, .11509E+00, .11689E+00, .11869E+00, .12047E+00,
          .12225E+00, .12403E+00, .12580E+00, .12758E+00, .12937E+00,
          .13117E+00, .13300E+00, .13486E+00, .13676E+00, .13871E+00,
          .14073E+00, .14284E+00, .14504E+00, .14736E+00, .14982E+00,
          .15242E+00, .15519E+00, .15813E+00, .16128E+00, .16464E+00,
          .16822E+00, .17205E+00, .17612E+00, .18045E+00, .18507E+00,
          .18999E+00, .19524E+00, .20085E+00, .20684E+00, .21324E+00,
          .22008E+00, .22741E+00, .23525E+00, .24364E+00, .25261E+00,
          .26219E+00, .27243E+00, .28339E+00, .29510E+00, .30760E+00,
          .32085E+00, .33506E+00, .35018E+00, .36608E+00, .38316E+00,
          .40125E+00, .42009E+00, .44023E+00, .46116E+00, .48364E+00,
          .50688E+00, .53174E+00, .55732E+00, .58462E+00, .61261E+00,
          .64172E+00, .67286E+00, .70447E+00, .73719E+00, .77223E+00,
          .80749E+00, .84384E+00, .88267E+00, .92153E+00, .96142E+00,
          .10023E+01, .10460E+01, .10893E+01, .11334E+01, .11784E+01,
          .12242E+01, .12730E+01, .13205E+01, .13686E+01, .14172E+01,
          .14661E+01, .15153E+01, .15672E+01, .16169E+01, .16664E+01,
          .17155E+01, .17642E+01, .18121E+01, .18591E+01, .19049E+01,
          .19495E+01, .19942E+01, .20355E+01, .20748E+01, .21118E+01,
          .21465E+01, .21786E+01, .22080E+01, .22345E+01, .22581E+01,
          .22787E+01, .22963E+01, .23109E+01, .23225E+01, .23311E+01,
          .23370E+01, .23401E+01, .23409E+01, .23393E+01, .23358E+01,
          .23305E+01, .23236E+01, .23156E+01, .23066E+01, .22969E+01,
          .22868E+01, .22766E+01, .22664E+01, .22564E+01, .22470E+01,
          .22382E+01, .22301E+01, .22230E+01, .22170E+01, .22120E+01,
          .22082E+01, .22056E+01, .22043E+01, .22042E+01, .22053E+01,
          .22077E+01, .22113E+01, .22162E+01, .22222E+01, .22294E+01,
          .22376E+01, .22470E+01, .22573E+01, .22687E+01, .22809E+01,
          .22941E+01, .23081E+01, .23229E+01, .23385E+01, .23548E+01,
          .23718E+01, .23894E+01, .24077E+01, .24265E+01, .24458E+01,
          .24656E+01, .24859E+01, .25066E+01, .25277E+01, .25493E+01,
          .25711E+01, .25934E+01, .26159E+01, .26387E+01
       };
    
      // linear interpolation to calculate kappa(T) W/(m.K)
      i_index = int(floor(fabs((T-T0)/dT)));
      Temp_Argon = T0+i_index*dT;
      kappa_T_Argon = kappa_Argon[i_index]
                    + (kappa_Argon[i_index+1]-kappa_Argon[i_index])*(T-Temp_Argon)/dT;

      return kappa_T_Argon;

   // end of kappa version implemented for argon plasma
}


template<class Thermo>
inline Foam::scalar Foam::arTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{
	//old version
    //return mu(p, T)*rPr_;
    //pr=nu/alpha=mu*cp/kappa
    //mu/pr = kappa/cp
    
    return kappa(p,T)/this->Cpv(p,T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::arTransport<Thermo>& Foam::arTransport<Thermo>::operator=
(
    const arTransport<Thermo>& ct
)
{
    Thermo::operator=(ct);

    mu_ = ct.mu_;
    rPr_ = ct.rPr_;

    return *this;
}


template<class Thermo>
inline void Foam::arTransport<Thermo>::operator+=
(
    const arTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(st);

    if (mag(molr1) + mag(st.nMoles()) > SMALL)
    {
        molr1 /= this->nMoles();
        scalar molr2 = st.nMoles()/this->nMoles();

        mu_ = molr1*mu_ + molr2*st.mu_;
        rPr_ = 1.0/(molr1/rPr_ + molr2/st.rPr_);
    }
}


template<class Thermo>
inline void Foam::arTransport<Thermo>::operator-=
(
    const arTransport<Thermo>& st
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(st);

    if (mag(molr1) + mag(st.nMoles()) > SMALL)
    {
        molr1 /= this->nMoles();
        scalar molr2 = st.nMoles()/this->nMoles();

        mu_ = molr1*mu_ - molr2*st.mu_;
        rPr_ = 1.0/(molr1/rPr_ - molr2/st.rPr_);
    }
}


template<class Thermo>
inline void Foam::arTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::arTransport<Thermo> Foam::operator+
(
    const arTransport<Thermo>& ct1,
    const arTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) + static_cast<const Thermo&>(ct2)
    );

    if (mag(ct1.nMoles()) + mag(ct2.nMoles()) < SMALL)
    {
        return arTransport<Thermo>
        (
            t,
            0,
            ct1.rPr_
        );
    }
    else
    {
        scalar molr1 = ct1.nMoles()/t.nMoles();
        scalar molr2 = ct2.nMoles()/t.nMoles();

        return arTransport<Thermo>
        (
            t,
            molr1*ct1.mu_ + molr2*ct2.mu_,
            1.0/(molr1/ct1.rPr_ + molr2/ct2.rPr_)
        );
    }
}


template<class Thermo>
inline Foam::arTransport<Thermo> Foam::operator-
(
    const arTransport<Thermo>& ct1,
    const arTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) - static_cast<const Thermo&>(ct2)
    );

    if (mag(ct1.nMoles()) + mag(ct2.nMoles()) < SMALL)
    {
        return arTransport<Thermo>
        (
            t,
            0,
            ct1.rPr_
        );
    }
    else
    {
        scalar molr1 = ct1.nMoles()/t.nMoles();
        scalar molr2 = ct2.nMoles()/t.nMoles();

        return arTransport<Thermo>
        (
            t,
            molr1*ct1.mu_ - molr2*ct2.mu_,
            1.0/(molr1/ct1.rPr_ - molr2/ct2.rPr_)
        );
    }
}


template<class Thermo>
inline Foam::arTransport<Thermo> Foam::operator*
(
    const scalar s,
    const arTransport<Thermo>& ct
)
{
    return arTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
        ct.mu_,
        1.0/ct.rPr_
    );
}


template<class Thermo>
inline Foam::arTransport<Thermo> Foam::operator==
(
    const arTransport<Thermo>& ct1,
    const arTransport<Thermo>& ct2
)
{
    return ct2 - ct1;
}


// ************************************************************************* //
