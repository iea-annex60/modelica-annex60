within IBPSA.Media.Refrigerants.Interfaces.Choices;
type InputChoice = enumeration(
    pT "(p,T) as inputs",
    dT "(d,T) as inputs",
    ph "(p,h) as inputs",
    ps "(p,s) as inputs")
   "Enumeration to define input choice of calculating thermodynamic state"
    annotation (Evaluate=true);
