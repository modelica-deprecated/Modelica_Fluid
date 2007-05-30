package Junctions "Junction components" 
  extends Modelica_Fluid.Icons.VariantLibrary;
  
  model JunctionIdeal 
    "Splitting/joining component with static balances for an infinitesimal control volume" 
    import Modelica_Fluid.Types;
    import Modelica_Fluid.Types.PortFlowDirection;
    
    replaceable package Medium=Modelica.Media.Interfaces.PartialMedium 
      "Fluid medium model" 
      annotation (choicesAllMatching=true);
    
    Modelica_Fluid.Interfaces.FluidPort_b port_1(
      redeclare package Medium=Medium,
      m_flow(min=if (portFlowDirection_1==PortFlowDirection.Entering) then 0.0 else -Modelica.Constants.inf,
      max=if (portFlowDirection_1==PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) 
      annotation (extent=[-120,-10; -100,10]);
    Modelica_Fluid.Interfaces.FluidPort_b port_2(
      redeclare package Medium=Medium,
      m_flow(min=if (portFlowDirection_2==PortFlowDirection.Entering) then 0.0 else -Modelica.Constants.inf,
      max=if (portFlowDirection_2==PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) 
      annotation (extent=[100,-10; 120,10]);
    Modelica_Fluid.Interfaces.FluidPort_b port_3(
      redeclare package Medium=Medium,
      m_flow(min=if (portFlowDirection_3==PortFlowDirection.Entering) then 0.0 else -Modelica.Constants.inf,
      max=if (portFlowDirection_3==PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) 
      annotation (extent=[-10,100; 10,120]);
    
    Medium.AbsolutePressure p "Pressure";
    Medium.SpecificEnthalpy h "Mixing enthalpy";
    Medium.MassFraction Xi[Medium.nXi] 
      "Independent mixture mass fractions m_i/m";
    Medium.ExtraProperty C[Medium.nC] "Trace substance mixture content";
    
    parameter Types.Init.Temp initType=Types.Init.NoInit 
      "Initialization option" 
      annotation(Evaluate=true,Dialog(tab="Initialization"));
    parameter Medium.AbsolutePressure p_start "Start value of pressure" 
      annotation(Dialog(tab="Initialization"));
    parameter Boolean use_T_start=true "=true, use T_start, otherwise h_start" 
      annotation(Dialog(tab="Initialization"),Evaluate=true);
    parameter Medium.Temperature T_start=
      if use_T_start then Medium.T_default else Medium.temperature_phX(p_start,h_start,X_start) 
      "Start value of temperature" 
      annotation(Dialog(tab="Initialization",enable=use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.specificEnthalpy_pTX(p_start,T_start,X_start) else Medium.h_default 
      "Start value of specific enthalpy" 
      annotation(Dialog(tab="Initialization",enable=not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX]=Medium.X_default 
      "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization",enable=Medium.nXi>0));
    
    parameter PortFlowDirection.Temp portFlowDirection_1=PortFlowDirection.Bidirectional 
      "Flow direction for port_1" 
     annotation(Dialog(tab="Advanced"));
    parameter PortFlowDirection.Temp portFlowDirection_2=PortFlowDirection.Bidirectional 
      "Flow direction for port_2" 
     annotation(Dialog(tab="Advanced"));
    parameter PortFlowDirection.Temp portFlowDirection_3=PortFlowDirection.Bidirectional 
      "Flow direction for port_3" 
     annotation(Dialog(tab="Advanced"));
    
    annotation(Documentation(info="<html>
  This model is the simplest implementation for a splitting/joining component for
  three flows. Its use is not required. It just formulates the balance
  equations in the same way that the connect symmantics would formulate them anyways.
  The main advantage of using this component is, that the user does not get
  confused when looking at the specific enthalpy at each port which might be confusing
  when not using a splitting/joining component. The reason for the confusion is that one exmanins the mixing
  enthalpy of the infinitesimal control volume introduced with the connect statement when
  looking at the specific enthalpy in the connector which
  might not be equal to the specific enthalpy at the port in the \"real world\".</html>"),
      Coordsys(grid=[1,1]),
      Icon(
        Rectangle(extent=[-100,41; 100,-47],   style(
            color=0,
            gradient=2,
            fillColor=8)),
        Rectangle(extent=[-100,37; 100,-43],   style(
            color=69,
            gradient=2,
            fillColor=69)),
        Rectangle(extent=[-34,100; 34,37], style(
            color=0,
            rgbcolor={0,0,0},
            gradient=1,
            fillColor=8,
            rgbfillColor={192,192,192})),
        Rectangle(extent=[-30,100; 30,35], style(
            color=69,
            rgbcolor={0,127,255},
            gradient=1,
            fillColor=69,
            rgbfillColor={0,127,255}))));
  initial equation 
    // Initial conditions
    if initType == Types.Init.NoInit then
      // no initial equations
    elseif initType == Types.Init.InitialValues then
      p = p_start;
      h = h_start;
    elseif initType == Types.Init.SteadyState then
      der(p) = 0;
      der(h) = 0;
    elseif initType == Types.Init.SteadyStateHydraulic then
      der(p) = 0;
      h = h_start;
    else
      assert(false, "Unsupported initialization option");
    end if;
  equation 
    port_1.m_flow + port_2.m_flow + port_3.m_flow = 0 "Mass balance";
    port_1.mXi_flow + port_2.mXi_flow + port_3.mXi_flow = zeros(Medium.nXi) 
      "Component mass balances";
    port_1.mC_flow + port_2.mC_flow + port_3.mC_flow = zeros(Medium.nC) 
      "Trace substance mass balances";
    
    port_1.H_flow = semiLinear(port_1.m_flow, port_1.h, h);
    port_2.H_flow = semiLinear(port_2.m_flow, port_2.h, h);
    port_3.H_flow = semiLinear(port_3.m_flow, port_3.h, h);
    
    port_1.mXi_flow = semiLinear(port_1.m_flow, port_1.Xi, Xi);
    port_2.mXi_flow = semiLinear(port_2.m_flow, port_2.Xi, Xi);
    port_3.mXi_flow = semiLinear(port_3.m_flow, port_3.Xi, Xi);
    
    port_1.mC_flow = semiLinear(port_1.m_flow, port_1.C, C);
    port_2.mC_flow = semiLinear(port_2.m_flow, port_2.C, C);
    port_3.mC_flow = semiLinear(port_3.m_flow, port_3.C, C);
    
    // Momentum balances
    port_1.p = p;
    port_2.p = p;
    port_3.p = p;
    
    port_1.H_flow + port_2.H_flow + port_3.H_flow = 0 "Energy balance";
  end JunctionIdeal;
  
  annotation (Documentation(info="<html>
 
</html>"));
  model JunctionVolume 
    "Splitting/joining component with static balances for a dynamic control volume" 
    import Modelica_Fluid.Types;
    import Modelica_Fluid.Types.PortFlowDirection;
    
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
      "Fluid medium model" 
        annotation (choicesAllMatching=true);
    parameter SI.Volume V "Volume";
    
    SI.InternalEnergy U "Internal energy";
    SI.Mass m "Total mass";
    SI.Mass[Medium.nXi] mXi "Independent masses";
    
    Modelica_Fluid.Interfaces.FluidPort_b port_1(
      redeclare package Medium=Medium,
      m_flow(min=if (portFlowDirection_1==PortFlowDirection.Entering) then 0.0 else -Modelica.Constants.inf,
      max=if (portFlowDirection_1==PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) 
      annotation (extent=[-120,-10; -100,10]);
    Modelica_Fluid.Interfaces.FluidPort_b port_2(
      redeclare package Medium=Medium,
      m_flow(min=if (portFlowDirection_2==PortFlowDirection.Entering) then 0.0 else -Modelica.Constants.inf,
      max=if (portFlowDirection_2==PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) 
      annotation (extent=[100,-10; 120,10]);
    Modelica_Fluid.Interfaces.FluidPort_b port_3(
      redeclare package Medium=Medium,
      m_flow(min=if (portFlowDirection_3==PortFlowDirection.Entering) then 0.0 else -Modelica.Constants.inf,
      max=if (portFlowDirection_3==PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) 
      annotation (extent=[-10,100; 10,120]);
    
    Medium.ExtraProperty C[Medium.nC] "Trace substance mixture content";
    Medium.BaseProperties medium;
    
    parameter Types.Init.Temp initType=Types.Init.NoInit 
      "Initialization option" 
      annotation(Evaluate=true,Dialog(tab="Initialization"));
    parameter Medium.AbsolutePressure p_start "Start value of pressure" 
      annotation(Dialog(tab="Initialization"));
    parameter Boolean use_T_start=true "=true, use T_start, otherwise h_start" 
      annotation(Dialog(tab="Initialization"),Evaluate=true);
    parameter Medium.Temperature T_start=
      if use_T_start then Medium.T_default else Medium.temperature_phX(p_start,h_start,X_start) 
      "Start value of temperature" 
      annotation(Dialog(tab="Initialization",enable=use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.specificEnthalpy_pTX(p_start,T_start,X_start) else Medium.h_default 
      "Start value of specific enthalpy" 
      annotation(Dialog(tab="Initialization",enable=not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX]=Medium.X_default 
      "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization",enable=Medium.nXi>0));
    
    parameter PortFlowDirection.Temp portFlowDirection_1=PortFlowDirection.Bidirectional 
      "Flow direction for port_1" 
     annotation(Dialog(tab="Advanced"));
    parameter PortFlowDirection.Temp portFlowDirection_2=PortFlowDirection.Bidirectional 
      "Flow direction for port_2" 
     annotation(Dialog(tab="Advanced"));
    parameter PortFlowDirection.Temp portFlowDirection_3=PortFlowDirection.Bidirectional 
      "Flow direction for port_3" 
     annotation(Dialog(tab="Advanced"));
    
    annotation (Icon(
        Rectangle(extent=[-100,41; 100,-47],   style(
            color=0,
            gradient=2,
            fillColor=8)),
        Rectangle(extent=[-100,37; 100,-43],   style(
            color=69,
            gradient=2,
            fillColor=69)),
        Rectangle(extent=[-34,100; 34,37], style(
            color=0,
            rgbcolor={0,0,0},
            gradient=1,
            fillColor=8,
            rgbfillColor={192,192,192})),
        Rectangle(extent=[-30,100; 30,35], style(
            color=69,
            rgbcolor={0,127,255},
            gradient=1,
            fillColor=69,
            rgbfillColor={0,127,255})),
        Ellipse(extent=[-9,10; 11,-10],   style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=0,
              rgbfillColor={0,0,0}))),
      Coordsys(grid=[1,1]),
      Diagram);
  initial equation 
    // Initial conditions
    if initType == Types.Init.NoInit then
      // no initial equations
    elseif initType == Types.Init.InitialValues then
      medium.p = p_start;
      medium.h = h_start;
    elseif initType == Types.Init.SteadyState then
      der(medium.p) = 0;
      der(medium.h) = 0;
    elseif initType == Types.Init.SteadyStateHydraulic then
      der(medium.p) = 0;
      medium.h = h_start;
    else
      assert(false, "Unsupported initialization option");
    end if;
    
  equation 
    port_1.m_flow + port_2.m_flow + port_3.m_flow = der(m) "Mass balance";
    port_1.mXi_flow + port_2.mXi_flow + port_3.mXi_flow = der(mXi) 
      "Component mass balances";
    port_1.mC_flow + port_2.mC_flow + port_3.mC_flow = zeros(Medium.nC) 
      "Trace substance mass balances";
    
    port_1.H_flow = semiLinear(port_1.m_flow, port_1.h,medium.h);
    port_2.H_flow = semiLinear(port_2.m_flow, port_2.h, medium.h);
    port_3.H_flow = semiLinear(port_3.m_flow, port_3.h, medium.h);
    
    port_1.mXi_flow = semiLinear(port_1.m_flow, port_1.Xi, medium.Xi);
    port_2.mXi_flow = semiLinear(port_2.m_flow, port_2.Xi, medium.Xi);
    port_3.mXi_flow = semiLinear(port_3.m_flow, port_3.Xi, medium.Xi);
    
    port_1.mC_flow = semiLinear(port_1.m_flow, port_1.C, C);
    port_2.mC_flow = semiLinear(port_2.m_flow, port_2.C, C);
    port_3.mC_flow = semiLinear(port_3.m_flow, port_3.C, C);
    
    // Momentum balance
    //Suitable for compressible media
    port_1.p = medium.p;
    port_2.p = medium.p;
    port_3.p = medium.p;
    
    U=m*medium.u;
    mXi=m*medium.Xi;
    m=medium.d*V;
    
    port_1.H_flow + port_2.H_flow + port_3.H_flow = der(U) "Energy balance";
  end JunctionVolume;
  
  model MassFlowRatio "simple flow multiplier" 
    replaceable package Medium=Modelica.Media.Interfaces.PartialMedium annotation(choicesAllMatching);
    parameter Real ratio=1 "flow multiplier from port a to port b";
    Modelica_Fluid.Interfaces.FluidPort_a port_a(
                                  redeclare package Medium=Medium) 
      annotation (extent=[-110,-10; -90,10]);
    Modelica_Fluid.Interfaces.FluidPort_b port_b(
                                  redeclare package Medium = Medium) 
      annotation (extent=[90,-10; 110,10]);
  equation 
    port_b.m_flow = ratio*port_a.m_flow;
    port_b.H_flow = ratio*port_a.H_flow;
    port_b.mXi_flow = ratio*port_b.mXi_flow;
    port_b.mC_flow = ratio*port_b.mC_flow;
    
    port_a.p = port_b.p;
    port_a.H_flow = semiLinear(
      port_a.m_flow,
      port_a.h,
      port_b.h);
    port_a.mXi_flow = semiLinear(
      port_a.m_flow,
      port_a.Xi,
      port_b.Xi);
    port_a.mC_flow = semiLinear(
      port_a.m_flow,
      port_a.C,
      port_b.C);
    annotation (Icon(
        Line(points=[-80,0; 80,80], style(
            color=69,
            rgbcolor={0,128,255},
            thickness=4)),
        Line(points=[-80,0; 80,0], style(
            color=69,
            rgbcolor={0,128,255},
            thickness=4)),
        Line(points=[-80,0; 80,-80], style(
            color=69,
            rgbcolor={0,128,255},
            thickness=4)),
        Line(points=[-80,0; 80,40], style(
            color=69,
            rgbcolor={0,128,255},
            thickness=4)),
        Line(points=[-80,0; 80,-40], style(
            color=69,
            rgbcolor={0,128,255},
            thickness=4))), Documentation(info="<html>
This model reflects a simple flow multiplication, which is very helpful in cases where the flow is evenly distributed to several parallel flow paths which are identical in their dimensions and boundary conditions, as e.g. in heat exchangers. Only one of the parallel pipes needs to be simulated then. All flow variables in <b>port_b</b> are achieved by multiplying those in <b>port_a</b> with parameter <b>ratio</b>.
</html>"));
  end MassFlowRatio;
  
  model HeatFlowRatio "simple heat flow multiplier" 
    parameter Real ratio=1 "heat flow ratio from port_a to port_b";
    annotation (Icon(
        Line(points=[-80,0; 80,80], style(
            color=1,
            rgbcolor={255,0,0},
            thickness=4)),
        Line(points=[-80,0; 80,40], style(
            color=1,
            rgbcolor={255,0,0},
            thickness=4)),
        Line(points=[-80,0; 80,0], style(
            color=1,
            rgbcolor={255,0,0},
            thickness=4)),
        Line(points=[-80,0; 80,-40], style(
            color=1,
            rgbcolor={255,0,0},
            thickness=4)),
        Line(points=[-80,0; 80,-80], style(
            color=1,
            rgbcolor={255,0,0},
            thickness=4))), Documentation(info="<html>
Simple model for heat flow multiplication between the two ports. The heat flow rate in port_a is multiplied by parameter <b>ratio</b> to achieve the rate at port_b. Both temperatures are equal, the model may be used e.g. for parallel pipes in heat exchangers.
</html>"));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a port_a 
      annotation (extent=[-100,-10; -80,10]);
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b port_b 
      annotation (extent=[80,-10; 100,10]);
  equation 
    port_b.Q_flow=ratio*port_a.Q_flow;
    port_b.T=port_a.T;
  end HeatFlowRatio;
end Junctions;
