// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(U) 
  + fvm::div(phi, U)
  + turb.divDevRhoReff(U)
 ==
    fvOptions(U)
);

UEqn().relax();

fvOptions.constrain(UEqn());

if (momentumPredictor)
{
    solve(UEqn() == -fvc::grad(p));

    fvOptions.correct(U);
}
