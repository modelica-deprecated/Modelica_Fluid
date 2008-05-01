within Modelica_Fluid;
package Junctions "Junction components" 
  extends Modelica_Fluid.Icons.VariantLibrary;
  
  model JunctionIdeal 
    "Splitting/joining component with static balances for an infinitesimal control volume" 
    import Modelica_Fluid.Types;
    import Modelica_Fluid.Types.PortFlowDirection;
    
    replaceable package Medium=Modelica.Media.Interfaces.PartialMedium 
      "Fluid medium model" 
      annotation (choicesAllMatching=true);
    
    parameter PortFlowDirection.Temp portFlowDirection_1=PortFlowDirection.Bidirectional 
      "Flow direction for port_1" 
     annotation(Dialog(tab="Advanced"));
    parameter PortFlowDirection.Temp portFlowDirection_2=PortFlowDirection.Bidirectional 
      "Flow direction for port_2" 
     annotation(Dialog(tab="Advanced"));
    parameter PortFlowDirection.Temp portFlowDirection_3=PortFlowDirection.Bidirectional 
      "Flow direction for port_3" 
     annotation(Dialog(tab="Advanced"));
    
    Modelica_Fluid.Interfaces.FluidPort_a port_1(redeclare package Medium = 
          Medium, m_flow(min=if (portFlowDirection_1 == PortFlowDirection.Entering) then 
                  0.0 else -Modelica.Constants.inf, max=if (portFlowDirection_1
             == PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) 
      annotation (extent=[-110,-10; -90,10]);
    Modelica_Fluid.Interfaces.FluidPort_b port_2(redeclare package Medium = 
          Medium, m_flow(min=if (portFlowDirection_2 == PortFlowDirection.Entering) then 
                  0.0 else -Modelica.Constants.inf, max=if (portFlowDirection_2
             == PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) 
      annotation (extent=[90,-10; 110,10]);
    Modelica_Fluid.Interfaces.FluidPort_a port_3(
      redeclare package Medium=Medium,
      m_flow(min=if (portFlowDirection_3==PortFlowDirection.Entering) then 0.0 else -Modelica.Constants.inf,
      max=if (portFlowDirection_3==PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) 
      annotation (extent=[-10,90; 10,110]);
    
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
            rgbfillColor={0,127,255})),
        Text(
          extent=[-150,-60; 150,-100],
          style(color=3, rgbcolor={0,0,255}),
          string="%name")),
      Diagram);
    
  equation 
    connect(port_1, port_2) annotation (points=[-100,0; 100,0], style(
        color=69,
        rgbcolor={0,127,255},
        smooth=0));
    connect(port_1, port_3) annotation (points=[-100,0; 0,0; 0,100], style(
        color=69,
        rgbcolor={0,127,255},
        smooth=0));
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
    
    Modelica_Fluid.Interfaces.FluidPort_a port_1(
      redeclare package Medium=Medium,
      m_flow(min=if (portFlowDirection_1==PortFlowDirection.Entering) then 0.0 else -Modelica.Constants.inf,
      max=if (portFlowDirection_1==PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) 
      annotation (extent=[-110,-10; -90,10]);
    Modelica_Fluid.Interfaces.FluidPort_b port_2(
      redeclare package Medium=Medium,
      m_flow(min=if (portFlowDirection_2==PortFlowDirection.Entering) then 0.0 else -Modelica.Constants.inf,
      max=if (portFlowDirection_2==PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) 
      annotation (extent=[90,-10; 110,10]);
    Modelica_Fluid.Interfaces.FluidPort_a port_3(
      redeclare package Medium=Medium,
      m_flow(min=if (portFlowDirection_3==PortFlowDirection.Entering) then 0.0 else -Modelica.Constants.inf,
      max=if (portFlowDirection_3==PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf)) 
      annotation (extent=[-10,90; 10,110]);
    
    Medium.ExtraProperty C[Medium.nC] "Trace substance mixture content";
    Medium.BaseProperties medium(preferredMediumStates=true);
    
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
              rgbfillColor={0,0,0})),
        Text(
          extent=[-150,-60; 150,-100],
          style(color=3, rgbcolor={0,0,255}),
          string="%name")),
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
    // Boundary conditions
    port_1.h_outflow = medium.h;
    port_2.h_outflow = medium.h;
    port_3.h_outflow = medium.h;
    
    // Internal quantities
    m   = medium.d*V;
    mXi = m*medium.Xi;
    U   = m*medium.u;
    
    // Mass balances
    der(m)   = port_1.m_flow + port_2.m_flow + port_3.m_flow "Mass balance";
    der(mXi) = streamFlow(port_1.Xi_outflow) +
               streamFlow(port_2.Xi_outflow) +
               streamFlow(port_3.Xi_outflow) "Component mass balances";
    
    zeros(Medium.nC) = streamFlow(port_1.C_outflow) +
                       streamFlow(port_2.C_outflow) +
                       streamFlow(port_3.C_outflow) 
      "Trace substance mass balances";
    
    // Momentum balance (suitable for compressible media)
    port_1.p = medium.p;
    port_2.p = medium.p;
    port_3.p = medium.p;
    
    // Energy balance
    der(U) = streamFlow(port_1.h_outflow) +
             streamFlow(port_2.h_outflow) +
             streamFlow(port_3.h_outflow);
  end JunctionVolume;
  
  model MassFlowRatio "simple flow multiplier" 
    replaceable package Medium=Modelica.Media.Interfaces.PartialMedium annotation(choicesAllMatching);
    parameter Integer nOutlets=1 
      "Number of outlet ports (mass is distributed evenly between the outlet ports";
    Modelica_Fluid.Interfaces.FluidPort_a port_a(
                                  redeclare package Medium=Medium) 
      annotation (extent=[-110,-10; -90,10]);
    Modelica_Fluid.Interfaces.FluidPorts_b ports_b[nOutlets] 
                                    annotation (extent=[90,-40; 110,40]);
    
    annotation (Icon(
        Line(points=[-80,0; 80,0], style(
            color=69,
            rgbcolor={0,128,255},
            thickness=4)),
        Line(points=[-80,0; 80,28], style(
            color=69,
            rgbcolor={0,128,255},
            thickness=4)),
        Line(points=[-80,0; 80,-28], style(
            color=69,
            rgbcolor={0,128,255},
            thickness=4)),
        Text(
          extent=[-150,100; 150,60],
          style(color=3, rgbcolor={0,0,255}),
          string="%name")), Documentation(info="<html>
<p>
This model describes a simple flow partitioning, which is very helpful in cases where the flow is evenly distributed to several parallel flow paths which are identical in their dimensions and boundary conditions, as e.g. in heat exchangers. Only one of the parallel pipes needs to be simulated then. All flow variables in <b>port_b[i]</b> are equal to those at <b>port_a</b> divided by <b>nOutlets</b>. All effort variables are equal at all ports.
</p>
</html>"),
      Diagram);
    
  equation 
    port_a.h_outflow  =  sum( {inflow(ports_b[i].h_outflow)     for i in 1:nOutlets})/nOutlets;
    port_a.Xi_outflow = {sum( {inflow(ports_b[i].Xi_outflow[j]) for i in 1:nOutlets})/nOutlets for j in 1:Medium.nXi};
    port_a.C_outflow  = {sum( {inflow(ports_b[i].C_outflow[j])  for i in 1:nOutlets})/nOutlets for j in 1:Medium.nXi};
    
    for i in 1:nOutlets loop
       ports_b[i].h_outflow  = inflow(port_a.h_outflow);
       ports_b[i].Xi_outflow = inflow(port_a.Xi_outflow);
       ports_b[i].C_outflow  = inflow(port_a.C_outflow);
      
       // Momentum balance  
       port_a.p = ports_b[i].p;
      
       // Mass balance
       ports_b[i].m_flow = -port_a.m_flow/nOutlets;
    end for;
  end MassFlowRatio;
  
  model HeatFlowRatio "simple heat flow multiplier" 
    parameter Integer nOutlets=1 
      "Number of outlet ports (heat is distributed evenly between the outlet ports";
    annotation (Icon(
        Line(points=[-80,0; 80,30], style(
            color=42,
            rgbcolor={127,0,0},
            thickness=4)),
        Line(points=[-80,0; 80,0], style(
            color=42,
            rgbcolor={127,0,0},
            thickness=4)),
        Line(points=[-80,0; 80,-32], style(
            color=42,
            rgbcolor={127,0,0},
            thickness=4)),
        Text(
          extent=[-150,100; 150,60],
          style(color=3, rgbcolor={0,0,255}),
          string="%name")), Documentation(info="<html>
<p>
Simple model for heat flow partitioning between the two ports. The heat flow rate in port_a is divided by parameter <b>nOutlets</b> to achieve the rate at ports port_b. All temperatures are equal. The model may be used e.g. for parallel pipes in heat exchangers.
</p>
</html>"));
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort_a 
      annotation (extent=[-110,-10; -90,10]);
    Modelica_Fluid.Interfaces.HeatPorts_b heatPorts_b[nOutlets] 
                                   annotation (extent=[90,-40; 110,40]);
  equation 
    for i in 1:nOutlets loop
       heatPorts_b[i].Q_flow = -heatPort_a.Q_flow/nOutlets;
       heatPorts_b[i].T      =  heatPort_a.T;
    end for;
  end HeatFlowRatio;
  
  model GenericJunction 
    "Branching component with balances for a dynamic control volume" 
    import Modelica_Fluid.Types;
    import Modelica_Fluid.Types.PortFlowDirection;
    import Modelica_Fluid.Types.ModelStructure;
    parameter Integer n_a=1 "Number of ports on side a";
    parameter Integer n_b=1 "Number of ports on side b";
    replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
      "Fluid medium model" 
        annotation (choicesAllMatching=true);
    parameter SI.Volume V "Volume";
    parameter SI.Pressure dp_nom "nominal (linear) pressure drop" annotation(Dialog(enable=not modelStructure==ModelStructure.avb));
    parameter SI.MassFlowRate mflow_nom "nominal mass flow rate"  annotation(Dialog(enable=not modelStructure==ModelStructure.avb));
    
    SI.InternalEnergy U "Internal energy";
    SI.Mass m "Total mass";
    SI.Mass[Medium.nXi] mXi "Independent masses";
    
    Interfaces.FluidStatePorts_a[n_a] ports_a(
      redeclare each package Medium=Medium,
      m_flow(each min=if (portFlowDirection==PortFlowDirection.Entering) then 0.0 else -Modelica.Constants.inf,
      each max=if (portFlowDirection==PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf))  annotation (extent=[-110,40;
          -90,-40]);
    Interfaces.FluidStatePorts_b[n_b] ports_b(
      redeclare each package Medium=Medium,
      m_flow(each min=if (portFlowDirection==PortFlowDirection.Entering) then 0.0 else -Modelica.Constants.inf,
      each max=if (portFlowDirection==PortFlowDirection.Leaving) then 0.0 else Modelica.Constants.inf))  annotation (extent=[90,40;
          110,-40]);
    Medium.ExtraProperty C[Medium.nC] "Trace substance mixture content";
    Medium.BaseProperties medium(T(start=T_start),p(start=p_start),h(start=h_start),X(start=X_start), preferredMediumStates=true);
    
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
    
    parameter PortFlowDirection.Temp portFlowDirection=PortFlowDirection.Bidirectional 
      "Allowed flow direction for ports" 
     annotation(Dialog(tab="Advanced"));
    parameter ModelStructure.Temp modelStructure=ModelStructure.avb annotation(Evaluate=true);
    
    Medium.EnthalpyFlowRate ports_a_H_flow[n_a];
    Medium.EnthalpyFlowRate ports_b_H_flow[n_b];
    Medium.MassFlowRate ports_a_mXi_flow[n_a,Medium.nXi];
    Medium.MassFlowRate ports_b_mXi_flow[n_b,Medium.nXi];
    Medium.ExtraPropertyFlowRate ports_a_mC_flow[n_a,Medium.nC];
    Medium.ExtraPropertyFlowRate ports_b_mC_flow[n_b,Medium.nC];
    
    annotation (Icon(
        Ellipse(extent=[-19,0; 1,-20],    style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=0,
              rgbfillColor={0,0,0})),
        Ellipse(extent=[-100,100; 100,-100],
                                         style(
            color=10,
            rgbcolor={135,135,135},
            gradient=3,
            fillColor=69,
            rgbfillColor={0,128,255})),
        Ellipse(extent=[-9,10; 11,-10],   style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=0,
              rgbfillColor={0,0,0})),Text(
        extent=[-150,150; 150,110],
        string="%name",
        style(gradient=2, fillColor=69))),
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
    sum(ports_a.m_flow)+sum(ports_b.m_flow)=der(m) "Mass balance";
    
    sum(ports_a_H_flow) + sum(ports_b_H_flow) = der(U) "Energy balance";
    
    for i in 1:Medium.nXi loop
      sum(ports_a_mXi_flow[i,:])+sum(ports_b_mXi_flow[i,:]) = der(mXi[i]);
    end for;
    
    for i in 1:Medium.nC loop
      sum(ports_a_mC_flow[i,:])+sum(ports_b_mC_flow[i,:]) = 0;
    end for;
    
    for i in 1:n_a loop
      ports_a[i].h_outflow  = medium.h;
      ports_a[i].Xi_outflow = medium.Xi;
      ports_a[i].C_outflow = C;
      
      ports_a_H_flow[i] = semiLinear(ports_a[i].m_flow, inflow(ports_a[i].h_outflow), medium.h) 
        "Energy balance";
      ports_a_mXi_flow[i,:] = semiLinear(ports_a[i].m_flow, inflow(ports_a[i].Xi_outflow), medium.Xi) 
        "Component mass balances";
      ports_a_mC_flow[i,:] = semiLinear(ports_a[i].m_flow, inflow(ports_a[i].C_outflow), C) 
        "Trace substance mass balances";
    end for;
    
    for i in 1:n_b loop
      ports_b[i].h_outflow  = medium.h;
      ports_b[i].Xi_outflow = medium.Xi;
      ports_b[i].C_outflow = C;
      
      ports_b_H_flow[i] = semiLinear(ports_b[i].m_flow, inflow(ports_b[i].h_outflow), medium.h) 
        "Energy balance";
      ports_b_mXi_flow[i,:] = semiLinear(ports_b[i].m_flow, inflow(ports_b[i].Xi_outflow), medium.Xi) 
        "Component mass balances";
      ports_b_mC_flow[i,:] = semiLinear(ports_b[i].m_flow, inflow(ports_b[i].C_outflow), C) 
        "Trace substance mass balances";
    end for;
    
    if modelStructure==ModelStructure.avb or modelStructure == ModelStructure.av_b then
      ports_a.p=fill(medium.p, n_a);
    else
      ports_a.p-fill(medium.p,n_a) = ports_a.m_flow*dp_nom/mflow_nom;
    end if;
    
    if modelStructure==ModelStructure.avb or modelStructure==ModelStructure.a_vb then
      ports_b.p=fill(medium.p,n_b);
    else
      ports_b.p-fill(medium.p,n_b)=ports_b.m_flow*dp_nom/mflow_nom;
    end if;
    
    U=m*medium.u;
    mXi=m*medium.Xi;
    m=medium.d*V;
    
  end GenericJunction;
end Junctions;
