package Obsolete "Components that are obsolete and will be removed" 
  extends Modelica.Icons.Library;
  package Types 
    package InitTypes 
      "Type, constants and menu choices to define initialization, as temporary solution until enumerations are available" 
      extends Modelica.Icons.Enumeration;
      constant Integer NoInit = 0 "No explicit initial conditions";
      constant Integer InitialValues = 1 
        "Initial conditions specified by initial values";
      constant Integer SteadyState = 2 "Full steady-state";
      constant Integer SteadyStateHydraulic = 3 
        "Hydraulic steady state, initial values of other state variables given";
      type Temp 
        "Temporary type with choices for menus (until enumerations are available)" 
        extends Integer(min=0,max=3);
        annotation (Evaluate=true, choices(
            choice=Modelica_Fluid.Types.Init.NoInit 
              "NoInit (No explicit initial conditions)",
            choice=Modelica_Fluid.Types.Init.InitialValues 
              "InitialValues (Initial conditions specified by initial values)",
            choice=Modelica_Fluid.Types.Init.SteadyState 
              "SteadyState (Full steady state)",
            choice=Modelica_Fluid.Types.Init.SteadyStateHydraulic 
              "SteadyStateHydraulic (Hydraulic steady state, initial values of other states given)"));
      end Temp;
    end InitTypes;
  end Types;
  
end Obsolete;
