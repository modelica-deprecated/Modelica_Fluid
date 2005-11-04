package Types 
  package Init 
    "Type, constants and menu choices to define initialization, as temporary solution until enumerations are available" 
    annotation (preferedView="text");
    extends Modelica.Icons.Enumeration;
    constant Integer NoInit = 1 
      "No initial conditions (guess values for p, T or h, X)";
    constant Integer InitialValues = 2 "Initial values for p, T or h, X";
    constant Integer SteadyState = 3 
      "Steady state (guess values for p, T or h, X)";
    constant Integer SteadyStateHydraulic = 4 
      "Hydraulic steady state (der(p)=0), guess value for p, initial values for T or h, X";
    type Temp 
      "Temporary type with choices for menus (until enumerations are available)" 
      extends Integer(min=1, max=4);
      annotation (Evaluate=true, choices(
          choice=Modelica_Fluid.WorkInProgress.Types.Init.NoInit 
            "NoInit (guess values for p, T or h, X)",
          choice=Modelica_Fluid.WorkInProgress.Types.Init.InitialValues 
            "InitialValues (initial values for p, T or h, X)",
          choice=Modelica_Fluid.WorkInProgress.Types.Init.SteadyState 
            "SteadyState (guess values for p, T or h, X)",
          choice=Modelica_Fluid.WorkInProgress.Types.Init.SteadyStateHydraulic 
            "SteadyStateHydraulic (der(p)=0, guess value for p, initial values for T or h, X)"));
    end Temp;
  end Init;

  package InitWithGlobalDefault 
    "Type, constants and menu choices to define initialization, as temporary solution until enumerations are available" 
    annotation (preferedView="text");
    extends Modelica.Icons.Enumeration;
    constant Integer NoInit = 1 
      "No initial conditions (guess values for p, T or h, X)";
    constant Integer InitialValues = 2 "Initial values for p, T or h, X";
    constant Integer SteadyState = 3 
      "Steady state (guess values for p, T or h, X)";
    constant Integer SteadyStateHydraulic = 4 
      "Hydraulic steady state (der(p)=0), guess value for p, initial values for T or h, X";
    constant Integer UseEnvironmentOption = 5 
      "Use initialization defined in environment component";
    
    type Temp 
      "Temporary type with choices for menus (until enumerations are available)" 
      extends Integer(min=1, max=5);
      annotation (Evaluate=true, choices(
          choice=Modelica_Fluid.WorkInProgress.Types.Init.NoInit 
            "NoInit (guess values for p, T or h, X)",
          choice=Modelica_Fluid.WorkInProgress.Types.Init.InitialValues 
            "InitialValues (initial values for p, T or h, X)",
          choice=Modelica_Fluid.WorkInProgress.Types.Init.SteadyState 
            "SteadyState (guess values for p, T or h, X)",
          choice=Modelica_Fluid.WorkInProgress.Types.Init.SteadyStateHydraulic 
            "SteadyStateHydraulic (der(p)=0, guess value for p, initial values for T or h, X)",
          choice=Modelica_Fluid.WorkInProgress.Types.Init.UseEnvironmentOption 
            "UseEnvironmentOption (use initialization defined in environment component)"));
    end Temp;
  end InitWithGlobalDefault;

  package Flow 
    "Type, constants and menu choices to define whether flow reversal is allowed, as temporary solution until enumerations are available" 
    annotation (preferedView="text");
    extends Modelica.Icons.Enumeration;
    constant Integer Unidirectional = 1 "Fluid flows only in one direction";
    constant Integer Bidirectional = 2 
      "No restrictions on fluid flow (flow reversal possible)";
    
    type Temp 
      "Temporary type with choices for menus (until enumerations are available)" 
      extends Integer(min=1, max=2);
      annotation (Evaluate=true, choices(
          choice=Modelica_Fluid.WorkInProgress.Types.FlowReversal.Unidirectional 
            "Unidirectional (fluid flows only in one direction)",
          choice=Modelica_Fluid.WorkInProgress.Types.FlowReversal.Bidirectional 
            "Bidirectional (flow reversal possible)"));
    end Temp;
  end Flow;

  package FlowWithGlobalDefault 
    "Type, constants and menu choices to define whether flow reversal is allowed, as temporary solution until enumerations are available" 
    annotation (preferedView="text");
    extends Modelica.Icons.Enumeration;
    constant Integer Unidirectional = 1 "Fluid flows only in one direction";
    constant Integer Bidirectional = 2 
      "No restrictions on fluid flow (flow reversal possible)";
    constant Integer UseEnvironmentOption = 3 
      "Use FlowReversal defined in environment component";
    
    type Temp 
      "Temporary type with choices for menus (until enumerations are available)" 
      extends Integer(min=1, max=3);
      annotation (Evaluate=true, choices(
          choice=Modelica_Fluid.WorkInProgress.Types.FlowReversal.Unidirectional 
            "Unidirectional (fluid flows only in one direction)",
          choice=Modelica_Fluid.WorkInProgress.Types.FlowReversal.Bidirectional 
            "Bidirectional (flow reversal possible)",
          choice=Modelica_Fluid.WorkInProgress.Types.FlowReversal.UseEnvironmentOption 
            "UseEnvironmentOption (use FlowReversal defined in environment component)"));
    end Temp;
  end FlowWithGlobalDefault;
end Types;
