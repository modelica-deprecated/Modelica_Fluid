within ;
package Modelica_Fluid_WorkInProgress 
  "Contains models under development and related examples" 
  extends Modelica.Icons.Library;
  
annotation (Documentation(info="<html>
 
</html>"), uses(
      Modelica(version="2.2.2"),
      Modelica_Fluid(version="1.0 Stream Connectors Beta 1"),
      UserInteraction(version="0.53")));
  package Components 
    
  model IsolatedPipe 
      "Model of an isolated pipe consisting of n pipe segments/FiniteVolumes" 
      
    replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching = true);
      
    extends 
        Modelica_Fluid_WorkInProgress.Interfaces.PartialInitializationParameters;
      
    parameter Integer nVolumes(min=1)=1 
        "Number of pipe segments/finite volumes";
      
    parameter Modelica.SIunits.Length L "Length of pipe";
    parameter Modelica.SIunits.AbsolutePressure dp_nominal(min=1.e-10)=1 
        "|frictionType = ConstantLaminar or ConstantTurbulent| Nominal pressure drop";
      
    parameter Modelica.SIunits.MassFlowRate m_flow_nominal=1E-3 
        "Nominal mass flow rate at nominal pressure drop";
      
    parameter Modelica.SIunits.Area A_a;
    parameter Modelica.SIunits.Area A_b=A_a;
      
    parameter Modelica.SIunits.Length Z_a=0;
    parameter Modelica.SIunits.Length Z_b=Z_a;
      
    parameter Boolean dynamicMomentumBalance=false 
        "If false, der(m_flow) is neglected in momentum balance" 
                                                   annotation(Evaluate=true,
        Dialog(tab="Level of Detail"));
    parameter Boolean includeKineticTerm=false 
        "If false, d*v^2 is neglected in momentum balance" 
                                               annotation(Evaluate=true,
        Dialog(tab="Level of Detail"));
    parameter Boolean includeViscosity=false 
        "If false, artifical viscosity is neglected" 
                                            annotation(Evaluate=true, Dialog(tab=
            "Level of Detail"));
    parameter Real viscosityFactor1=0 annotation(Dialog(enable=includeViscosity,tab="Level of Detail"));
    parameter Real viscosityFactor2=1 annotation(Dialog(enable=includeViscosity,tab="Level of Detail"));
      
    Modelica_Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = Medium) 
                annotation (extent=[-120, -10; -100, 10]);
    Modelica_Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = Medium) 
                annotation (extent=[120, -10; 100, 10]);
    Modelica_Fluid_WorkInProgress.Utilities.PipeSegment pipeSegment[nVolumes](
        redeclare package Medium = Medium,
        each initType=initType,
        each p_start=p_start,
        each use_T_start=use_T_start,
        each T_start=T_start,
        each h_start=h_start,
        each X_start=X_start,
        each L=L/nVolumes,
        each dp_nominal=dp_nominal/nVolumes,
        each A_a=A_a "has to be corrected: linear distribution of A",
        each Z_a=Z_a "has to be corrected: linear distribution of Z",
        each m_flow_nominal=m_flow_nominal,
        each dynamicMomentumBalance=dynamicMomentumBalance,
        each includeKineticTerm=includeKineticTerm,
        each includeViscosity=includeViscosity,
        each viscosityFactor1=viscosityFactor1,
        each viscosityFactor2=viscosityFactor2);
      
  annotation (Icon(
      Rectangle(extent=[-100, 60; 100, -60], style(color=0, fillColor=8)),
      Rectangle(extent=[-100, 34; 100, -36], style(
          color=69,
          gradient=2,
          fillColor=69)),
      Text(
        extent=[-150, 125; 150, 65],
        string="%name",
        style(gradient=2, fillColor=69)),
        Ellipse(extent=[-58,14; -28,-14], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=0,
              rgbfillColor={0,0,0})),
        Ellipse(extent=[22,14; 52,-14], style(
              color=0,
              rgbcolor={0,0,0},
              fillColor=0,
              rgbfillColor={0,0,0}))));
  equation 
    connect(port_a, pipeSegment[1].port_a);
    connect(port_b, pipeSegment[nVolumes].port_b);
    for i in 1:nVolumes - 1 loop
      connect(pipeSegment[i].port_b, pipeSegment[i + 1].port_a);
    end for;
  end IsolatedPipe;
    
  model OpenTank "Tank with three inlet/outlet-arrays at variable heights" 
      import Modelica_Fluid.Types;
      import Modelica_Fluid;
    replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
      annotation (choicesAllMatching=true);
      
    parameter Modelica.SIunits.Height height "Height of tank";
    parameter Modelica.SIunits.Area area "Area of tank";
    parameter Modelica.SIunits.Volume V0=0 
        "Volume of the liquid when the level is zero";
    // parameter Real k = 4.9 "Heat transfer coefficient from tank to ambient";
      
    parameter Integer n_topPorts = 0 "Number of topPorts" 
      annotation(Dialog(group="topPorts (= pipes at top of tank; only flow into tank)"));
    parameter Modelica.SIunits.Height top_heights[n_topPorts]=fill(height,
          n_topPorts) "Heights of topPorts" 
      annotation(Dialog(group="topPorts (= pipes at top of tank; only flow into tank)",enable=n_topPorts > 0));
      
    parameter Integer n_bottomPorts = 0 "Number of bottomPorts" 
       annotation(Dialog(group="bottomPorts (= pipes at bottom of tank; in and out flow of tank)"));
    parameter Modelica.SIunits.Height bottom_heights[n_bottomPorts]=fill(0.0,
          n_bottomPorts) "Heights of bottomPorts" 
       annotation(Dialog(group="bottomPorts (= pipes at bottom of tank; in and out flow of tank)",enable=n_bottomPorts > 0));
    parameter Modelica.SIunits.Diameter bottom_diameters[n_bottomPorts] 
        "Inner (hydraulic) diameters of bottomPorts" 
       annotation(Dialog(group="bottomPorts (= pipes at bottom of tank; in and out flow of tank)", enable=n_bottomPorts > 0));
      
    parameter Integer n_sidePorts = 0 "Number of sidePorts" 
       annotation(Dialog(group="sidePorts (= pipes at side of tank; in and out flow of tank)"));
    parameter Modelica.SIunits.Height side_heights[n_sidePorts] 
        "Heights of sidePorts" 
       annotation(Dialog(group="sidePorts (= pipes at side of tank; in and out flow of tank)",enable=n_sidePorts > 0));
    parameter Modelica.SIunits.Area side_diameters[n_sidePorts] 
        "Inner (hydraulic) diameters of sidePorts" 
       annotation(Dialog(group="sidePorts (= pipes at side of tank; in and out flow of tank)",enable=n_sidePorts > 0));
      
    parameter Medium.AbsolutePressure p_ambient=ambient.default_p_ambient 
        "Tank surface pressure" 
      annotation(Dialog(tab = "Ambient and Initialization", group = "Ambient"));
    parameter Medium.Temperature T_ambient = ambient.default_T_ambient 
        "Tank surface Temperature" 
      annotation(Dialog(tab = "Ambient and Initialization", group = "Ambient"));
      
    parameter Types.Init.Temp initType=
              Types.Init.InitialValues "Initialization option" 
      annotation(Evaluate=true,Dialog(tab = "Ambient and Initialization", group = "Initialization"));
    parameter Modelica.SIunits.Height level_start "Start value of tank level" 
      annotation(Dialog(tab="Ambient and Initialization", group = "Initialization"));
    parameter Medium.Temperature T_start=T_ambient "Start value of temperature"
      annotation(Dialog(tab = "Ambient and Initialization", group = "Initialization"));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
        "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Ambient and Initialization", group = "Initialization", enable=Medium.nXi > 0));
      
  /*
  parameter Real k_small(min=0) = 0 
    "Small regularization range if tank level is below bottom_heights[i] or side_heights[i]; k_small = 0 gives ideal switch"
              annotation(Evaluate=true, Dialog(tab="Advanced"));
*/
      
    Modelica_Fluid.WorkInProgress.FluidStorage.FluidPort_ArrayIcon topPorts[
                                                                         n_topPorts](redeclare 
          package Medium=Medium, m_flow(each start=0), mXi_flow(each start=0)) 
      annotation (extent=[-30,100; 30,108]);
    Modelica_Fluid.WorkInProgress.FluidStorage.FluidPort_ArrayIcon bottomPorts[
                                                                            n_bottomPorts](redeclare 
          package Medium=Medium, m_flow(each start=0), mXi_flow(each start=0)) 
      annotation (extent=[-30,-108; 30,-100],   rotation=90);
    Modelica_Fluid.WorkInProgress.FluidStorage.FluidPort_ArrayIcon sidePorts[
                                                                          n_sidePorts](redeclare 
          package Medium=Medium, m_flow(each start=0), mXi_flow(each start=0)) 
      annotation (extent=[100,30; 108,-30]);
      
    Medium.BaseProperties medium(
      preferredMediumStates=true,
      p(start=p_ambient),
      T(start=T_start),
      h(start=h_start),
      Xi(start=X_start[1:Medium.nXi]));
      
    Modelica.SIunits.Height level(stateSelect=StateSelect.prefer, start=
            level_start) "Level height of tank";
    Modelica.SIunits.Volume V(stateSelect=StateSelect.never) 
        "Actual tank volume";
    Modelica.SIunits.Energy U "Internal energy of tank volume";
    Modelica.SIunits.Mass m "Mass of fluid in tank";
    Modelica.SIunits.Mass mXi[Medium.nXi] 
        "Masses of independent components in the fluid";
      
  Modelica_Fluid.Ambient ambient;
    protected 
    Real Q_lost = 0 "Wärmeverlust (currently dummy)";
    parameter Medium.SpecificEnthalpy h_start = Medium.specificEnthalpy_pTX(p_ambient, T_start, X_start);
    parameter Integer precision = 3 "Precision for tank level in animation" annotation(Hide=false);
      
    Medium.EnthalpyFlowRate H_flow_topPorts[n_topPorts];
    Medium.EnthalpyFlowRate port_b_H_flowottomPorts[n_bottomPorts];
    Medium.EnthalpyFlowRate H_flow_sidePorts[n_sidePorts];
      
    Medium.MassFlowRate m_flow_topPorts[n_topPorts];
    Medium.MassFlowRate m_flow_bottomPorts[n_bottomPorts];
    Medium.MassFlowRate m_flow_sidePorts[n_sidePorts];
      
    Medium.MassFlowRate mXi_flow_topPorts[n_topPorts,Medium.nXi];
    Medium.MassFlowRate port_b_mXi_flowottomPorts[n_bottomPorts,Medium.nXi];
    Medium.MassFlowRate mXi_flow_sidePorts[n_sidePorts,Medium.nXi];
      
    protected 
    Modelica_Fluid.WorkInProgress.FluidStorage.TankAttachment tankAttachmentTop[
                                                                             n_topPorts](
      each h=medium.h,
      each d=medium.d,
      each Xi=medium.Xi,
      each p_ambient=p_ambient,
      each level=level,
      each h_start = h_start,
      each X_start = X_start,
      each level_start = level_start,
      each onlyInFlow=true,
      H_flow=H_flow_topPorts,
      m_flow=m_flow_topPorts,
      mXi_flow=mXi_flow_topPorts,
      pipeHeight=top_heights,
      redeclare package Medium = Medium) 
        annotation (extent=[-20,80; 20,40]);
    Modelica_Fluid.WorkInProgress.FluidStorage.TankAttachment 
        tankAttachmentBottom[                                                   n_bottomPorts](
      each h=medium.h,
      each d=medium.d,
      each Xi=medium.Xi,
      each p_ambient=p_ambient,
      each level=level,
      each h_start = h_start,
      each X_start = X_start,
      each level_start = level_start,
      each onlyInFlow=false,
      H_flow=port_b_H_flowottomPorts,
      m_flow=m_flow_bottomPorts,
      mXi_flow=port_b_mXi_flowottomPorts,
      pipeHeight=bottom_heights,
      pipeDiameter=bottom_diameters,
      redeclare package Medium = Medium) 
        annotation (extent=[-20,-80; 20,-40]);
    Modelica_Fluid.WorkInProgress.FluidStorage.TankAttachment 
        tankAttachmentSide[                                                   n_sidePorts](
      each h=medium.h,
      each d=medium.d,
      each Xi=medium.Xi,
      each p_ambient=p_ambient,
      each level=level,
      each h_start = h_start,
      each X_start = X_start,
      each level_start = level_start,
      each onlyInFlow=false,
      H_flow=H_flow_sidePorts,
      m_flow=m_flow_sidePorts,
      mXi_flow=mXi_flow_sidePorts,
      pipeHeight=side_heights,
      pipeDiameter=side_diameters,
      redeclare package Medium = Medium) 
        annotation (extent=[40,-20; 80,20], rotation=90);
      
  equation 
    for i in 1:n_bottomPorts loop
      connect(tankAttachmentBottom[i].port, bottomPorts[i]) 
        annotation (points=[0,-80.4; 0,-104],
                                            style(color=3, rgbcolor={0,0,255}));
    end for;
      
    for i in 1:n_topPorts loop
       connect(tankAttachmentTop[i].port, topPorts[i]) 
         annotation (points=[0,80.4; 0,104],
                                           style(color=3, rgbcolor={0,0,255}));
    end for;
      
    for i in 1:n_sidePorts loop
       connect(tankAttachmentSide[i].port, sidePorts[i]) 
         annotation (points=[80.4,-1.24914e-015; 104,0],
                                             style(color=3, rgbcolor={0,0,255}));
    end for;
      
    medium.p = p_ambient;
    V = area*level+V0 "Volume of fluid";
    m = V*medium.d "Mass of fluid";
    mXi = m*medium.Xi "Mass of fluid components";
    U = m*medium.u "Internal energy of fluid";
    // Q_lost = - k*2*sqrt(Modelica.Constants.pi*area)*level*(medium.T - T_ambient);
      
    // Mass balances
    der(m) = sum(m_flow_bottomPorts) + sum(m_flow_sidePorts) + sum(m_flow_topPorts);
    for i in 1:Medium.nXi loop
         der(mXi[i]) = sum(port_b_mXi_flowottomPorts[i,:]) +
                       sum(mXi_flow_sidePorts[i,:]) +
                       sum(mXi_flow_topPorts[i,:]);
    end for;
      
    // Energy balance
    if Medium.singleState then
       der(U) = sum(port_b_H_flowottomPorts)+sum(H_flow_sidePorts)+sum(H_flow_topPorts) + Q_lost 
          "Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)";
    else
       der(U) = sum(port_b_H_flowottomPorts)+sum(H_flow_sidePorts)+sum(H_flow_topPorts) - p_ambient*der(V) + Q_lost;
    end if;
    assert(level <= height, "Tank is full (level = height = " + String(level) + ")");
      
  initial equation 
    if initType == Types.Init.NoInit then
      // no initial equations
    elseif initType == Types.Init.InitialValues then
      level = level_start;
      medium.T = T_start;
      medium.Xi = X_start[1:Medium.nXi];
    elseif initType == Types.Init.SteadyState then
      der(level) = 0;
      der(medium.h) = 0;
      der(medium.Xi) = zeros(Medium.nXi);
    elseif initType == Types.Init.SteadyStateHydraulic then
      der(level) = 0;
      medium.T = T_start;
      medium.Xi = X_start[1:Medium.nXi];
    else
      assert(false, "Unsupported initialization option");
    end if;
    annotation (defaultComponentName="tank",
      Icon(
        Rectangle(extent=[-100,100; 100,0],     style(color=7, fillColor=7)),
        Rectangle(extent=DynamicSelect([-100,-100; 100,10],
                                       [-100, -100; 100, (-100 + 200*level/height)]),
            style(
            color=69,
            rgbcolor={0,127,255},
            fillColor=71,
            rgbfillColor={85,170,255},
            fillPattern=1)),
        Line(points=[-100, 100; -100, -100; 100, -100; 100, 100], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=69,
            rgbfillColor={0,127,255},
            fillPattern=1)),
        Text(
          extent=[-95,90; 95,60],
          string="%name",
          style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=69,
            rgbfillColor={0,127,255},
            fillPattern=1)),
        Text(
          extent=[-95,50; 95,35],
          style(color=0),
          string="%level_start",
        Line(points=[-100,100; 100,100], style(
            color=0,
            rgbcolor={0,0,0},
            pattern=3))),
        Text(
          extent=[-95,30; 95,5],
          style(color=0),
          string=DynamicSelect(" ",realString(level,1,integer(precision)))),
        Line(points=[-100,100; 100,100], style(
            color=0,
            rgbcolor={0,0,0},
            pattern=3)),
        Text(
          extent=[-40,124; -22,112],
          string="1",
          style(color=10, rgbcolor={135,135,135})),
        Text(
          extent=[111,37; 129,25],
          string="1",
          style(color=10, rgbcolor={135,135,135})),
        Text(
          extent=[-38,-112; -20,-124],
          string="1",
          style(color=10, rgbcolor={135,135,135}))),
      Documentation(info="<HTML>
<p>
This is a simplified model of a tank. 
The top part is open to the environment at the fixed pressure 
<tt>p_ambient</tt>. Heat transfer to the environment and to 
the tank walls is neglected.
The tank is filled with a single or multiple-substance liquid, 
assumed to have uniform temperature and mass fractions.
</p>
 
<p>
There are three arrays of FluidConnectors (topPorts, bottomPorts, sidePorts)
at which other fluid components can be attached. The difference between
these ports is just the graphical layout on the icon, to visualize
whether a pipe is connected to the top, bottom or side of the tank.
All ports are defined by height (with respect to bottom of tank)
and the diameter of the component that is attached to this port.
</p>
 
<p>
In a diagram animation, the fill level of the tank
is visualized, as well as the value of level.
</p>
 
</HTML>", revisions="<html>
<ul>
<li><i>Jan. 6, 2006</i>
   implemented by Katja Poschlad, Manuel Remelhel (AST Uni Dortmund), 
   Martin Otter (DLR), based on Components.Tank model.</li>
</ul>
</html>"),
      Diagram,
      uses(Modelica(version="2.2.1"), Modelica_Fluid(version="0.952")),
      Coordsys(grid=[1,1], scale=0.2));
  end OpenTank;
    
  model Tank "Obsolet component (use instead Components.OpenTank)" 
      import Modelica_Fluid.Types;
      import Modelica_Fluid;
    replaceable package Medium = 
        Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
      annotation (choicesAllMatching=true);
    parameter Modelica.SIunits.Area area "Tank area";
    parameter Modelica.SIunits.Area pipeArea "Area of outlet pipe";
    parameter Modelica.SIunits.Volume V0=0 
        "Volume of the liquid when the level is zero";
    parameter Modelica.SIunits.Height H0=0 
        "Height of zero level reference over the bottom port";
    parameter Medium.AbsolutePressure p_ambient=ambient.default_p_ambient 
        "Tank surface pressure";
    parameter Types.Init.Temp initType=
              Types.Init.NoInit "Initialization option" 
      annotation(Dialog(tab = "Initialization"));
    parameter Boolean use_T_start = true 
        "Use T_start if true, otherwise h_start" 
      annotation(Dialog(tab = "Initialization"), Evaluate = true);
    parameter Medium.Temperature T_start=
      if use_T_start then Medium.T_default else Medium.temperature_phX(p_ambient,h_start,X_start) 
        "Start value of temperature" 
      annotation(Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=
      if use_T_start then Medium.specificEnthalpy_pTX(p_ambient, T_start, X_start) else Medium.h_default 
        "Start value of specific enthalpy" 
      annotation(Dialog(tab = "Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
        "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    parameter Modelica.SIunits.Height level_start(min=0) "Initial tank level" 
      annotation(Dialog(tab="Initialization"));
      
    Modelica_Fluid.Interfaces.FluidPort_b port(
        redeclare package Medium = Medium,
        m_flow(start=0),
        mXi_flow(each start=0)) 
      annotation (extent=[-10, -120; 10, -100], rotation=90);
    Medium.BaseProperties medium(
      preferredMediumStates=true,
      p(start=p_ambient),
      T(start=T_start),
      Xi(start=X_start[1:Medium.nXi]));
      outer Modelica_Fluid.Ambient ambient "Ambient conditions";
      
    Modelica.SIunits.Height level(start=level_start, stateSelect=StateSelect.prefer) 
        "Height of tank level over the zero reference";
    Medium.AbsolutePressure p_bottom "Pressure at bottom of tank";
    Modelica.SIunits.Energy U "Internal energy of tank volume";
    Modelica.SIunits.Volume V(stateSelect=StateSelect.never) 
        "Actual tank volume";
    Real m(quantity=Medium.mediumName, unit="kg", stateSelect=StateSelect.never) 
        "Mass of tank volume";
    Real mX[Medium.nXi](quantity=Medium.substanceNames, each unit="kg") 
        "Component masses of the independent substances";
      
  equation 
    medium.p = p_ambient;
    V = area*level+V0 "Volume of fluid";
    m = V*medium.d "Mass of fluid";
    mX = m*medium.Xi "Mass of fluid components";
    U = m*medium.u "Internal energy of fluid";
      
    // Mass balance
    der(m) = port.m_flow;
    der(mX) = port.mXi_flow;
      
    // Momentum balance
    p_bottom = (medium.d*ambient.g*(level+H0)) + p_ambient;
    port.p = p_bottom - smooth(2,noEvent(if port.m_flow < 0 then port.m_flow^2/(2*medium.d*pipeArea^2) else 0));
      
    // Energy balance
    if Medium.singleState then
      der(U) = port.H_flow "Mechanical work is neglected";
    else
      der(U) = port.H_flow - p_ambient*der(V);
    end if;
      
    /* Handle reverse and zero flow */
    port.H_flow = semiLinear(port.m_flow, port.h, medium.h);
    port.mXi_flow = semiLinear(port.m_flow, port.Xi, medium.Xi);
      
  initial equation 
    if initType == Types.Init.NoInit then
      // no initial equations
    elseif initType == Types.Init.InitialValues then
      level = level_start;
      if use_T_start then
        medium.T = T_start;
      else
        medium.h = h_start;
      end if;
      medium.Xi = X_start[1:Medium.nXi];
    elseif initType == Types.Init.SteadyState then
      der(level) = 0;
      der(medium.h) = 0;
      der(medium.Xi) = zeros(Medium.nXi);
    elseif initType == Types.Init.SteadyStateHydraulic then
      der(level) = 0;
      if use_T_start then
        medium.T = T_start;
      else
        medium.h = h_start;
      end if;
      medium.Xi = X_start[1:Medium.nXi];
    else
      assert(false,"Unsupported initialization option initType = " + String(initType)
                   +"\nin model Modelica_Fluid.Components.Tank");
    end if;
    annotation (
      Icon(
        Rectangle(extent=[-100, 90; 100, 26], style(color=7, fillColor=7)),
        Rectangle(extent=[-100, 26; 100, -100], style(
            color=69,
            fillColor=69,
            fillPattern=1)),
        Line(points=[-100, 100; -100, -100; 100, -100; 100, 100], style(
            color=0,
            fillColor=69,
            fillPattern=1)),
        Text(
          extent=[-112, 162; 122, 102],
          string="%name",
          style(fillColor=69, fillPattern=1)),
        Text(
          extent=[-86, -38; 94, -78],
          style(color=0),
          string="%level_start"),
        Text(
          extent=[-94, 78; 94, 38],
          style(color=0),
          string="%p_ambient"),
        Text(
          extent=[-94, 14; 90, -2],
          style(color=0),
          string="level_start")),
      Documentation(info="<HTML>
<p>
This is a simplified model of a tank. The top part is open to the environment at the fixed pressure <tt>p_ambient</tt>. Heat transfer to the environment and to the tank walls is neglected.
The tank is filled with a single or multiple-substance liquid, assumed to have uniform temperature and mass fractions.
<p>The geometry of the tank is specified by the following parameters: <tt>V0</tt> is the volume of the liquid when the level is at the zero reference; <tt>area</tt> is the cross-sectional area of the tank; <tt>H0</tt> is the height of the zero-reference level plane over the port connector. It is thus possible to model rounded-bottom tanks, as long as they have a cylindrical shape in the range of operating levels.
<p>The tank can be initialized with the following options:
<ul>
<li>NoInit: no explicit initial conditions
<li>InitialValues: initial values of temperature (or specific enthalpy), composition and level are specified
<li>SteadyStateHydraulic: initial values of temperature (or specific enthalpy) and composition are specified; the initial level is determined so that levels and pressure are at steady state.
</ul>
Full steady state initialization is not supported, because the corresponding intial equations for temperature/enthalpy are undetermined (the flow rate through the port at steady state is zero). 
</p>
</HTML>",   revisions="<html>
<ul>
<li><i>1 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from a previous version of Modelica_Fluid</li>
</ul>
</html>"),
      Diagram);
  end Tank;
    
    model LumpedPipe 
      "Short pipe with one volume, wall friction and gravity effect" 
      replaceable package Medium = 
        Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching = true);
      parameter Modelica_Fluid.Types.Init.Temp initType=Modelica_Fluid.Types.Init.NoInit 
        "Initialization option" 
        annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p_a_start 
        "Start value of pressure at port_a" 
        annotation(Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p_b_start = p_a_start 
        "Start value of pressure at port_b" 
        annotation(Dialog(tab = "Initialization"));
      parameter Boolean use_T_start = true 
        "= true, use T_start, otherwise h_start" 
        annotation(Dialog(tab = "Initialization"), Evaluate=true);
      parameter Medium.Temperature T_start=
        if use_T_start then Medium.T_default else Medium.temperature_phX(p_a_start,h_start,X_start) 
        "Start value of temperature" 
        annotation(Dialog(tab = "Initialization", enable = use_T_start));
      parameter Medium.SpecificEnthalpy h_start=
        if use_T_start then Medium.specificEnthalpy_pTX(p_a_start, T_start, X_start) else Medium.h_default 
        "Start value of specific enthalpy" 
        annotation(Dialog(tab = "Initialization", enable = not use_T_start));
      parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
        "Start value of mass fractions m_i/m" 
        annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
      replaceable package WallFriction = 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.QuadraticTurbulent
        extends 
        Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.PartialWallFriction
        "Characteristic of wall friction"  annotation(choicesAllMatching=true);
      
      parameter Modelica.SIunits.Length length "Length of pipe";
      parameter Modelica.SIunits.Diameter diameter 
        "Inner (hydraulic) diameter of pipe";
      parameter Modelica.SIunits.Length height_ab=0.0 
        "Height(port_b) - Height(port_a)"                                annotation(Evaluate=true);
      parameter Modelica.SIunits.Length roughness(min=0)=2.5e-5 
        "Absolute roughness of pipe (default = smooth steel pipe)" 
          annotation(Dialog(enable=WallFriction.use_roughness));
      parameter Boolean use_nominal = false 
        "= true, if eta_nominal and d_nominal are used, otherwise computed from medium"
          annotation(Evaluate=true);
      parameter Modelica.SIunits.DynamicViscosity eta_nominal=0.01 
        "Nominal dynamic viscosity (for wall friction computation)" annotation(Dialog(enable=use_nominal));
      parameter Modelica.SIunits.Density d_nominal=0.01 
        "Nominal density (for wall friction computation)" annotation(Dialog(enable=use_nominal));
      parameter Modelica_Fluid.Types.FlowDirection.Temp flowDirection=
          Modelica_Fluid.Types.FlowDirection.Bidirectional 
        "Unidirectional (port_a -> port_b) or bidirectional flow component" 
         annotation(Dialog(tab="Advanced"));
      parameter Modelica.SIunits.AbsolutePressure dp_small=1 
        "Turbulent flow if |dp| >= dp_small (only used if WallFriction=QuadraticTurbulent)"
        annotation(Dialog(tab="Advanced", enable=WallFriction.use_dp_small));
      
        Modelica_Fluid.Interfaces.FluidPort_a port_a(
                                      redeclare package Medium = Medium) 
        "Fluid connector a (positive design flow direction is from port_a to port_b)"
          annotation (extent=[-110,-10; -90,10]);
        Modelica_Fluid.Interfaces.FluidPort_b port_b(
                                      redeclare package Medium = Medium) 
        "Fluid connector b (positive design flow direction is from port_a to port_b)"
          annotation (extent=[110,-10; 90,10]);
      
      annotation (defaultComponentName="pipe",Icon(
          Rectangle(extent=[-100,44; 100,-44],   style(
              color=0,
              gradient=2,
              fillColor=8)),
          Rectangle(extent=[-100,40; 100,-40],   style(
              color=69,
              gradient=2,
              fillColor=69)),
          Text(
            extent=[-145,-40; 155,-90],
            string="%name",
            style(gradient=2, fillColor=69)),
          Ellipse(extent=[-11,10; 9,-10],   style(
                color=0,
                rgbcolor={0,0,0},
                fillColor=0,
                rgbfillColor={0,0,0}))),       Documentation(info="<html>
<p>
Simple pipe model consisting of one volume, 
wall friction (with different friction correlations)
and gravity effect. This model is mostly used to demonstrate how
to build up more detailed models from the basic components.
Note, if the \"thermalPort\" is not connected, then the pipe
is totally insulated (= no thermal flow from the fluid to the
pipe wall/environment).
</p>
</html>"),
        Diagram,
        Coordsys(grid=[1,1], scale=0),
        uses(Modelica_Fluid(version="1.0 Beta 2"), Modelica(version="2.2.2")));
      Modelica_Fluid.PressureLosses.WallFrictionAndGravity frictionAndGravity1(
          redeclare package Medium = Medium,
          flowDirection=flowDirection,
          redeclare package WallFriction = WallFriction,
          length=length/2,
          diameter=diameter,
          height_ab=height_ab/2,
          roughness=roughness,
          use_nominal=use_nominal,
          eta_nominal=eta_nominal,
          d_nominal=d_nominal,
          from_dp=true,
          dp_small=dp_small,
          show_Re=false,
          p_a_start=p_a_start,
          p_b_start=(p_a_start+p_b_start)/2,
          T_start=T_start,
          h_start=h_start,
          X_start=X_start) annotation (extent=[-60,-10; -40,10]);
      Modelica_Fluid_WorkInProgress.Components.PortVolume volume(
        redeclare package Medium = Medium,
        V=Modelica.Constants.pi*(diameter/2)^2*length,
        initType=initType,
        p_start=(p_a_start + p_b_start)/2,
        use_T_start=use_T_start,
        T_start=T_start,
        h_start=h_start,
        X_start=X_start) 
        annotation (extent=[-10,-10; 10,10]);
      Modelica_Fluid.PressureLosses.WallFrictionAndGravity frictionAndGravity2(
          redeclare package Medium = Medium,
          flowDirection=flowDirection,
          redeclare package WallFriction = WallFriction,
          length=length/2,
          diameter=diameter,
          height_ab=height_ab/2,
          roughness=roughness,
          use_nominal=use_nominal,
          eta_nominal=eta_nominal,
          d_nominal=d_nominal,
          from_dp=true,
          dp_small=dp_small,
          show_Re=false,
          p_a_start=(p_a_start+p_b_start)/2,
          p_b_start=p_b_start,
          T_start=T_start,
          h_start=h_start,
          X_start=X_start,
          use_T_start=true) 
                           annotation (extent=[40,-10; 60,10]);
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a thermalPort 
        annotation (extent=[-10,44; 10,64]);
    equation 
      connect(frictionAndGravity1.port_a, port_a) 
        annotation (points=[-60,0; -100,0], style(color=69, rgbcolor={0,127,255}));
      connect(frictionAndGravity1.port_b, volume.port) 
        annotation (points=[-40,0; 0,0], style(color=69, rgbcolor={0,127,255}));
      connect(frictionAndGravity2.port_a, volume.port) 
        annotation (points=[40,0; 0,0], style(color=69, rgbcolor={0,127,255}));
      connect(frictionAndGravity2.port_b, port_b) 
        annotation (points=[60,0; 100,0], style(color=69, rgbcolor={0,127,255}));
      connect(volume.thermalPort, thermalPort) 
        annotation (points=[0,10; 0,54], style(color=42, rgbcolor={191,0,0}));
    end LumpedPipe;
    
      model PortVolume 
      "Obsolete (cannot be implemented with stream connectors) Fixed volume associated with a port by the finite volume method (used to build up physical components; fulfills mass and energy balance)" 
      import SI = Modelica.SIunits;
      import Modelica_Fluid.Types;
        annotation(__Dymola_obsolete="PortVolume cannot be implemented with stream connectors");
        extends 
        Modelica_Fluid_WorkInProgress.Interfaces.PartialInitializationParameters;
      
        replaceable package Medium = 
            Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
            annotation (choicesAllMatching = true);
      
        parameter SI.Volume V "Volume";
      
        Modelica_Fluid.Interfaces.FluidPort_a port(
          redeclare package Medium = Medium) "Fluid port" 
          annotation (extent=[-10, -10; 10, 10], rotation=0);
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a thermalPort 
        "Thermal port" 
          annotation (extent=[-10,90; 10,110]);
      
        Medium.BaseProperties medium(preferredMediumStates=true,
                    p(start=p_start), T(start=T_start),
                    h(start=h_start), Xi(start= X_start[1:Medium.nXi]));
        SI.Energy U "Internal energy of fluid";
        SI.Mass m "Mass of fluid";
        SI.Mass mXi[Medium.nXi] "Masses of independent components in the fluid";
        annotation (
         Icon(
            Ellipse(extent=[-100, 100; 100, -100], style(
                color=0,
                rgbcolor={0,0,0},
                gradient=3,
                fillColor=68,
                rgbfillColor={170,213,255})),
            Text(extent=[-150,-100; 150,-150], string="%name")),
                               Documentation(info="<html>
<p>
This component models the <b>volume</b> of <b>fixed size</b> that is
associated with the <b>fluid port</b> to which it is connected.
This means that all medium properties inside the volume, are identical
to the port medium properties. In particular, the specific enthalpy
inside the volume (= medium.h) is always identical to the specific enthalpy
in the port (port.h = medium.h). Usually, this model is used when
discretizing a component according to the finite volume method into
volumes in internal ports that only store energy and mass and into
transport elements that just transport energy, mass and momentum
between the internal ports without storing these quantities during the
transport. This splitting is only possible under certain assumptions.
</p>
</html>"),Diagram);
      equation 
        // medium properties set by port values
          port.p          = medium.p;
          port.h_outflow  = medium.h;
          port.Xi_outflow = medium.Xi;
          thermalPort.T   = medium.T;
      
        // Total quantities
           m    = V*medium.d "Total Mass";
           mXi = m*medium.Xi "Independent component masses";
           U    = m*medium.u "Internal energy";
      
        // Mass and energy balance
           der(m)    = port.m_flow "Total mass balance";
           der(mXi)  = semiLinear(port.m_flow,inflow(port.Xi_outflow),medium.Xi) 
        "Independent component mass balance";
           der(U)    = semiLinear(port.m_flow,inflow(port.h_outflow),medium.h) + thermalPort.Q_flow 
        "Energy balance";
           zeros(Medium.nC) = semiLinear(port.m_flow, inflow(port.C_outflow), port.C_outflow) 
        "Trace substances";
      
      initial equation 
        // Initial conditions
        if initType == Types.Init.NoInit then
          // no initial equations
        elseif initType == Types.Init.InitialValues then
          if not Medium.singleState then
             medium.p = p_start;
          end if;
          if use_T_start then
            medium.T = T_start;
          else
            medium.h = h_start;
          end if;
          medium.Xi = X_start[1:Medium.nXi];
        elseif initType == Types.Init.SteadyState then
          if not Medium.singleState then
             der(medium.p) = 0;
          end if;
          der(medium.h) = 0;
          der(medium.Xi) = zeros(Medium.nXi);
        elseif initType == Types.Init.SteadyStateHydraulic then
          if not Medium.singleState then
             der(medium.p) = 0;
          end if;
          if use_T_start then
            medium.T = T_start;
          else
            medium.h = h_start;
          end if;
          medium.Xi = X_start[1:Medium.nXi];
        else
          assert(false, "Unsupported initialization option");
        end if;
      end PortVolume;
  end Components;
  
  package Examples "Examples demonstrating the usage of components" 
    package ShockWaves "Examples demonstrating shock waves in pipes" 
    extends Modelica.Icons.Library;
      
    model TestSimpleIsolatedPipePressure "Test IsolatedPipe component" 
        
        import Modelica.SIunits.Conversions.*;
        
      extends Modelica.Icons.Example;
        
      parameter Real pressurePulseHeight=1E5;
      parameter Real pressurePulseWidth=50E-6;
      parameter Real pressurePulseStart=0.1E-3;
      parameter Real pressurePulseBase=2E5;
        
      Modelica_Fluid.Sources.FixedBoundary_pTX ambient(
        port(h(start=10000)),
          p=pressurePulseBase,
          T=320,
          redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
      annotation (extent=[80,50; 60,70]);
        
    annotation (
      Diagram,
      experiment(
            StopTime=0.0003,
            NumberOfIntervals=5000,
            fixedstepsize=1e-006,
            Algorithm=""),
      experimentSetupOutput,
      Commands(file="Simulate and plot pressure.mos", file=
         "Simulate and Plot Temperature.mos"));
        
      Modelica_Fluid_WorkInProgress.Components.IsolatedPipe isolatedPipe(
          redeclare package Medium = Modelica.Media.Air.DryAirNasa,
          dp_nominal=50,
          L=0.025,
          dynamicMomentumBalance=true,
          includeKineticTerm=true,
          nVolumes=25,
          p_start=2e5,
          T_start=320,
          A_a=1,
          includeViscosity=true,
          initType=Modelica_Fluid.Types.Init.InitialValues) 
               annotation (extent=[20,50; 40,70]);
      Modelica_Fluid.Sources.PrescribedBoundary_pTX prescribedAmbient(
                                                                     redeclare 
            package Medium=Modelica.Media.Air.DryAirNasa) 
                      annotation (extent=[-22,50; -2,70]);
      Modelica.Blocks.Math.Add Add1 
                                  annotation (extent=[-58,56; -38,76]);
      Modelica.Blocks.Sources.Ramp Ramp3(
          duration=scalar({pressurePulseWidth/2}),
          height=scalar({pressurePulseHeight}),
          offset=scalar({pressurePulseBase}),
          startTime=scalar({pressurePulseStart})) 
                     annotation (extent=[-100,30; -80,50]);
      Modelica.Blocks.Sources.Ramp Ramp4(
          duration=scalar({pressurePulseWidth/2}),
          height=scalar({-pressurePulseHeight}),
          offset=0,
          startTime=scalar({pressurePulseStart + pressurePulseWidth/2})) 
                      annotation (extent=[-100,90; -80,70]);
      Modelica.Blocks.Sources.Constant Constant1(k=320) 
        annotation (extent=[-60,20; -40,40]);
        
      UserInteraction.Outputs.SpatialPlot SpatialPlot1(
        y=isolatedPipe.pipeSegment.medium.p,
        x=linspace(0, 1, isolatedPipe.nVolumes),
          minY=pressurePulseBase - pressurePulseHeight,
          maxY=pressurePulseBase + pressurePulseHeight) 
                               annotation(extent=[-100,-100; 100,20]);
        
    equation 
      connect(isolatedPipe.port_b, ambient.port) 
      annotation (points=[41,60; 60,60], style(color=69));
      connect(prescribedAmbient.port, isolatedPipe.port_a)    annotation(points=[-2,60;
              19,60],      style(color=69, rgbcolor={0,127,255}));
      connect(Ramp3.y,Add1.u2)           annotation (points=[-79,40; -72,40; -72,60;
            -60,60],          style(color=3));
      connect(Ramp4.y,Add1.u1) 
      annotation (points=[-79,80; -70,80; -70,72; -60,72],   style(color=3));
        
      connect(Constant1.y, prescribedAmbient.T_in) annotation (points=[-39,30; -32,30;
            -32,60; -24,60], style(color=74, rgbcolor={0,0,127}));
        connect(Add1.y, prescribedAmbient.p_in) annotation (points=[-37,66; -24,66],
            style(color=74, rgbcolor={0,0,127}));
    end TestSimpleIsolatedPipePressure;
      
    model TestThreeIsolatedPipesPressure "Test ShortPipe componet" 
        
        import Modelica.SIunits.Conversions.*;
        
      extends Modelica.Icons.Example;
      parameter Real pressurePulsHeight=1E4;
      parameter Real pressurePulsWidth=1E-3;
      parameter Real pressurePulsStart=0.1E-3;
      parameter Modelica.Media.Interfaces.PartialMedium.Temperature T_ambient=300;
      parameter Modelica.Media.Interfaces.PartialMedium.Temperature T_ambient1=310;
        
    annotation (
      Diagram,
      experiment(
          StopTime=0.001,
          NumberOfIntervals=5000,
          Tolerance=1e-006,
          fixedstepsize=1e-006,
          Algorithm="Dassl"),
      experimentSetupOutput,
      Commands(file="Simulate and plot pressure.mos", file=
         "Simulate and Plot Temperature.mos"));
        
      Modelica_Fluid_WorkInProgress.Components.IsolatedPipe IsolatedPipe1(
          m_flow_nominal=1,
          A_a=0.01,
          dp_nominal=50,
          dynamicMomentumBalance=true,
          L=0.5,
          redeclare package Medium = Modelica.Media.Air.DryAirNasa,
          includeKineticTerm=true,
          includeViscosity=true,
          nVolumes=25,
          h_start=1e4,
          T_start=293.15,
          initType=Modelica_Fluid.Types.Init.InitialValues) 
               annotation (extent=[12,70; 32,90]);
      Modelica_Fluid.Sources.PrescribedBoundary_pTX prescribedAmbient(
          redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
                      annotation (extent=[-22,70; -2,90]);
      Modelica.Blocks.Sources.Ramp Ramp1(
          duration=scalar({pressurePulsWidth/2}),
          height=scalar({pressurePulsHeight}),
          offset=1E5,
          startTime=scalar({pressurePulsStart})) 
                     annotation (extent=[-100,52; -80,72]);
      Modelica.Blocks.Sources.Ramp Ramp2(
          duration=scalar({pressurePulsWidth/2}),
          height=scalar({-pressurePulsHeight}),
          offset=0,
          startTime=scalar({pressurePulsStart + pressurePulsWidth/2})) 
                      annotation (extent=[-100,106; -80,86]);
      Modelica.Blocks.Math.Add Add1 
                                  annotation (extent=[-60,76; -40,96]);
      Modelica_Fluid.Sources.FixedBoundary ambient1(
        p=1E5,
        port(h(start=10000)),
        T=T_ambient,
          redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
      annotation (extent=[100,70; 80,90]);
      Modelica_Fluid_WorkInProgress.Components.IsolatedPipe IsolatedPipe2(
          m_flow_nominal=1,
          dp_nominal=50,
          A_a=0.01/2,
          dynamicMomentumBalance=true,
          L=0.5,
          redeclare package Medium = Modelica.Media.Air.DryAirNasa,
          includeKineticTerm=true,
          includeViscosity=true,
          nVolumes=25,
          T_start=293.15,
          initType=Modelica_Fluid.Types.Init.InitialValues) 
               annotation (extent=[48,70; 68,90]);
      UserInteraction.Outputs.SpatialPlot SpatialPlot1(
        maxY=1.1e5,
        minY=0.9e5,
          x=linspace(0, 1, IsolatedPipe1.nVolumes),
        y=IsolatedPipe1.pipeSegment.medium.p) 
                            annotation(extent=[-100,-40; 0,20]);
      UserInteraction.Outputs.SpatialPlot SpatialPlot2(
        maxY=1.1e5,
        minY=0.9e5,
          x=linspace(0, 1, IsolatedPipe2.nVolumes),
        y=IsolatedPipe2.pipeSegment.medium.p) 
                            annotation(extent=[0,-40; 100,20]);
      Modelica_Fluid.Sources.FixedBoundary ambient2(
        p=1E5,
        port(h(start=10000)),
        T=T_ambient,
          redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
      annotation (extent=[100,30; 80,50]);
      Modelica_Fluid_WorkInProgress.Components.IsolatedPipe IsolatedPipe3(
          m_flow_nominal=1,
          dp_nominal=25,
          A_a=0.01/2,
          dynamicMomentumBalance=true,
          L=0.25,
          redeclare package Medium = Modelica.Media.Air.DryAirNasa,
          includeKineticTerm=true,
          includeViscosity=true,
          nVolumes=20,
          T_start=293.15,
          initType=Modelica_Fluid.Types.Init.InitialValues) 
               annotation (extent=[48,30; 68,50]);
      UserInteraction.Outputs.SpatialPlot SpatialPlot3(
        maxY=1.1e5,
        minY=0.9e5,
          x=linspace(0, 1, IsolatedPipe3.nVolumes),
        y=IsolatedPipe3.pipeSegment.medium.p) 
                            annotation(extent=[6,-100; 54,-40]);
      Modelica.Blocks.Sources.Constant Constant1(k=T_ambient1) 
        annotation (extent=[-60,46; -40,66]);
    equation 
        
      connect(prescribedAmbient.port, IsolatedPipe1.port_a) 
      annotation (points=[-2,80; 11,80], style(color=69));
      connect(Ramp1.y,Add1.u2)           annotation (points=[-79,62; -72,62; -72,80;
            -62,80],          style(color=3));
      connect(Ramp2.y,Add1.u1) 
      annotation (points=[-79,96; -70,96; -70,92; -62,92],   style(color=3));
      connect(IsolatedPipe2.port_a, IsolatedPipe1.port_b) 
                                                annotation(points=[47,80; 33,80],
          style(color=69, rgbcolor={0,127,255}));
      connect(IsolatedPipe2.port_b, ambient1.port) 
                                             annotation(points=[69,80; 80,80],
          style(color=69, rgbcolor={0,127,255}));
      connect(IsolatedPipe3.port_a, IsolatedPipe1.port_b) 
                                                annotation(points=[47,40; 40,40;
            40,80; 33,80], style(color=69, rgbcolor={0,127,255}));
      connect(IsolatedPipe3.port_b,ambient2. port) 
                                             annotation(points=[69,40; 80,40],
          style(color=69, rgbcolor={0,127,255}));
        
      connect(Add1.y, prescribedAmbient.p_in) annotation (points=[-39,86; -32,86;
            -32,86; -24,86], style(color=74, rgbcolor={0,0,127}));
      connect(Constant1.y, prescribedAmbient.T_in) annotation (points=[-39,56; -32,
            56; -32,80; -24,80], style(color=74, rgbcolor={0,0,127}));
    end TestThreeIsolatedPipesPressure;
      
    end ShockWaves;
    
    package Tanks "Examples with Tanks" 
      extends Modelica.Icons.Library;
      
      model TwoOpenTanks "Demonstrating the usage of OpenTank" 
        import Modelica_Fluid;
        extends Modelica.Icons.Example;
         // replaceable package Medium = Modelica.Media.Water.ConstantPropertyLiquidWater extends 
        // replaceable package Medium = Modelica.Media.Water.StandardWaterOnePhase extends 
        // replaceable package Medium = Modelica.Media.Incompressible.Examples.Glycol47 extends
         replaceable package Medium = 
            Modelica.Media.Water.ConstantPropertyLiquidWater                           extends 
          Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
            annotation (choicesAllMatching = true);
        
        Modelica_Fluid.WorkInProgress.Components.OpenTank tank1(
          height=10,
          area=1,
          n_bottomPorts=1,
          bottom_diameters={0.1},
          level_start=6,
          redeclare package Medium = Medium) 
                         annotation (extent=[-80,20; -40,60]);
        Modelica_Fluid.WorkInProgress.Components.OpenTank tank2(
          height=10,
          area=1,
          n_bottomPorts=1,
          bottom_diameters={0.1},
          level_start=3,
          redeclare package Medium = Medium) 
                         annotation (extent=[20,20; 60,60]);
        Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe1(
          length=1,
          p_a_start=ambient.default_p_ambient,
          p_b_start=ambient.default_p_ambient,
          T_start=ambient.default_T_ambient,
          diameter=0.1,
          height_ab=0,
          redeclare package Medium = Medium,
          redeclare package WallFriction = 
              Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.Detailed) 
                annotation (extent=[-20,-10; 0,10],  rotation=0);
        annotation (Diagram,
          experiment(StopTime=40),
          experimentSetupOutput,
          Documentation(info="<html>
  
</html>"));
        
        inner Modelica_Fluid.Ambient ambient 
                                         annotation (extent=[60,-34; 80,-14]);
      equation 
        connect(pipe1.port_a, tank1.bottomPorts[1]) 
                                                   annotation (points=[-20,0; -60,0;
              -60,19.2], style(color=69, rgbcolor={0,127,255}));
        connect(pipe1.port_b,tank2. bottomPorts[1]) annotation (points=[0,0; 40,0; 40,
              19.2], style(color=69, rgbcolor={0,127,255}));
      end TwoOpenTanks;
      
      model ThreeOpenTanks "Demonstrating the usage of OpenTank" 
        
        annotation (
          Diagram,
          experiment(StopTime=150),
          Coordsys(grid=[1, 1], component=[20, 20]),
          uses(Modelica_Fluid(version="0.952")),
          experimentSetupOutput,
          Documentation(info="<html> 
  
</html>"));
        
        extends Modelica.Icons.Example;
        // replaceable package Medium = Modelica.Media.Water.ConstantPropertyLiquidWater extends 
        // replaceable package Medium = Modelica.Media.Water.StandardWaterOnePhase extends 
        // replaceable package Medium = Modelica.Media.Incompressible.Examples.Glycol47 extends
         replaceable package Medium = 
           Modelica.Media.Water.ConstantPropertyLiquidWater                    extends 
          Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
            annotation (choicesAllMatching = true);
        
        Modelica_Fluid_WorkInProgress.Components.OpenTank tank1(
          area=1,
          V0=0,
          bottom_heights={0},
          redeclare package Medium = Medium,
          initType=Modelica_Fluid.Types.Init.InitialValues,
          height=10,
          bottom_diameters={0.05},
          n_sidePorts=2,
          side_heights={3,2},
          n_topPorts=0,
          n_bottomPorts=1,
          side_diameters={0.1,0.05},
          level_start=1) 
                        annotation (extent=[-40,0; 0,40]);
        
        Modelica_Fluid_WorkInProgress.Components.OpenTank tank2(
          area=1,
          V0=0,
          bottom_heights={0},
          bottom_diameters={0.1},
          redeclare package Medium = Medium,
          initType=Modelica_Fluid.Types.Init.InitialValues,
          height=10,
          level_start=9,
          n_bottomPorts=1) 
                         annotation (extent=[30,60; 70,100]);
        
        Modelica_Fluid_WorkInProgress.Components.OpenTank tank3(
          area=1,
          V0=0,
          level_start=6,
          top_heights={10},
          redeclare package Medium = Medium,
          initType=Modelica_Fluid.Types.Init.InitialValues,
          n_sidePorts=1,
          side_diameters={0.05},
          n_topPorts=1,
          side_heights={6.5},
          height=20)                                            annotation (extent=[-40,-90;
              0,-50]);
        
        Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe1(
          redeclare package Medium = Medium,
          length=1,
          height_ab=2,
          diameter=0.05,
          p_a_start=ambient.default_p_ambient,
          p_b_start=ambient.default_p_ambient,
          T_start=ambient.default_T_ambient,
          redeclare package WallFriction = 
              Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.NoFriction)
                annotation (extent=[-30,-40; -10,-20],
                                                     rotation=90);
        Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe2(
          redeclare package Medium = Medium,
          length=1,
          diameter=0.1,
          height_ab=2,
          p_a_start=ambient.default_p_ambient,
          p_b_start=ambient.default_p_ambient,
          T_start=ambient.default_T_ambient,
          redeclare package WallFriction = 
              Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.NoFriction)
                annotation (extent=[40,30; 60,50],   rotation=90);
        inner Modelica_Fluid.Ambient ambient 
          annotation (extent=[-90,-90; -70,-70]);
        Modelica_Fluid.PressureLosses.WallFrictionAndGravity pipe3(
          redeclare package Medium = Medium,
          length=1,
          height_ab=2,
          diameter=0.05,
          p_a_start=ambient.default_p_ambient,
          p_b_start=ambient.default_p_ambient,
          T_start=ambient.default_T_ambient,
          redeclare package WallFriction = 
              Modelica_Fluid.PressureLosses.BaseClasses.WallFriction.NoFriction)
                annotation (extent=[20,-40; 40,-20], rotation=90);
      equation 
        connect(tank2.bottomPorts[1], pipe2.port_b) 
                                          annotation (points=[50,59.2; 50,50],
            style(color=3, rgbcolor={0,0,255}));
        connect(tank1.bottomPorts[1], pipe1.port_b) 
                                          annotation (points=[-20,-0.8; -20,-20],
            style(color=3, rgbcolor={0,0,255}));
        connect(pipe1.port_a, tank3.topPorts[1]) 
                           annotation (points=[-20,-40; -20,-49.2],
            style(color=69, rgbcolor={0,127,255}));
        connect(pipe2.port_a, tank1.sidePorts[1]) annotation (points=[50,30; 50,23;
              0.8,23], style(color=69, rgbcolor={0,127,255}));
        connect(pipe3.port_b, tank1.sidePorts[2]) annotation (points=[30,-20; 30,16;
              0.8,16; 0.8,17], style(color=69, rgbcolor={0,127,255}));
        connect(pipe3.port_a, tank3.sidePorts[1]) annotation (points=[30,-40; 30,
              -70; 0.8,-70],
                        style(color=69, rgbcolor={0,127,255}));
      end ThreeOpenTanks;
    end Tanks;
    
    model TestIsolatedPipe "Test ShortPipe component" 
      import Modelica.SIunits.Conversions.*;
      
      extends Modelica.Icons.Example;
      annotation (
        Diagram(Text(
            extent=[-80,44; -34,36],
            string="water",
            style(color=0, rgbcolor={0,0,0})), Text(
            extent=[-82,-49; -36,-57],
            style(color=0, rgbcolor={0,0,0}),
            string="DetailedAir")),
        experiment(StopTime=3),
        Coordsys(grid=[1, 1], component=[20, 20]));
      Modelica_Fluid_WorkInProgress.Components.IsolatedPipe insulatedPipe_1(
        L=1,
        dp_nominal=1.e4,
        A_a=0.02*0.02,
        A_b=0.02*0.02,
        T_start=Modelica.SIunits.Conversions.from_degC(20),
        m_flow_nominal=0.1,
        nVolumes=10,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        initType=Modelica_Fluid.Types.Init.InitialValues) 
        annotation (extent=[-9,0; 11,20]);
      
      Modelica_Fluid.Sources.PrescribedMassFlowRate_TX massFlowSource_1(T=
            from_degC(30), redeclare package Medium = 
            Modelica.Media.Air.DryAirNasa) 
        annotation (extent=[-50,0; -30,20]);
      Modelica_Fluid.Sources.FixedBoundary_pTX ambient_1(
                                                        T=from_degC(15),
          redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
        annotation (extent=[50,0; 30,20]);
      Modelica.Blocks.Sources.Ramp massFlowSignal(
        height=3,
        offset=0.01,
        duration=3)  annotation (extent=[-90,0; -70,20]);
      Modelica_Fluid_WorkInProgress.Components.IsolatedPipe insulatedPipe_2(
        L=1,
        dp_nominal=1.e4,
        A_a=0.02*0.02,
        A_b=0.02*0.02,
        T_start=Modelica.SIunits.Conversions.from_degC(20),
        m_flow_nominal=0.1,
        nVolumes=10,
        redeclare package Medium = Modelica.Media.Air.DryAirNasa,
        initType=Modelica_Fluid.Types.Init.InitialValues) 
        annotation (extent=[-9,-40; 11,-20]);
      Modelica_Fluid.Sources.PrescribedMassFlowRate_TX massFlowSource_2(T=
            from_degC(30), redeclare package Medium = 
            Modelica.Media.Air.DryAirNasa) 
        annotation (extent=[-50,-40; -30,-20]);
      Modelica_Fluid.Sources.FixedBoundary_pTX ambient_2(
                                                        T=from_degC(15),
          redeclare package Medium = Modelica.Media.Air.DryAirNasa) 
        annotation (extent=[50,-40; 30,-20]);
    equation 
      connect(massFlowSource_1.port, insulatedPipe_1.port_a) 
        annotation (points=[-30,10; -10,10],  style(color=69));
      connect(insulatedPipe_1.port_b, ambient_1.port) 
        annotation (points=[12,10; 30,10],   style(color=69));
      connect(massFlowSource_2.port, insulatedPipe_2.port_a) 
        annotation (points=[-30,-30; -10,-30],style(color=69));
      connect(insulatedPipe_2.port_b, ambient_2.port) 
        annotation (points=[12,-30; 30,-30], style(color=69));
      connect(massFlowSignal.y, massFlowSource_1.m_flow_in) annotation (points=[-69,
            10; -63,10; -63,16; -49.3,16], style(color=74, rgbcolor={0,0,127}));
      connect(massFlowSignal.y, massFlowSource_2.m_flow_in) annotation (points=[-69,
            10; -63,10; -63,-24; -49.3,-24], style(color=74, rgbcolor={0,0,127}));
    end TestIsolatedPipe;
    
    model TestPressureLoss 
      extends Modelica.Icons.Example;
      replaceable package Medium = 
          Modelica.Media.Water.ConstantPropertyLiquidWater 
        extends Modelica.Media.Interfaces.PartialMedium 
        "Medium in all components"                        annotation (
        choicesAllMatching =                                                                            true);
      
      parameter 
        Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.LossFactorData
        lossFactors1(
          zeta1=0.5,
          D_a=0.1,
          D_b=0.1,
          Re_turbulent=100,
          D_Re=0.1,
        zeta2=1) "Loss factor without data for laminar region"   annotation (extent=[-100,-40; -80,-20]);
      parameter 
        Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.LossFactorData
        lossFactors2(
          zeta1=0.5,
          zetaLaminarKnown=true,
          c0=200,
          D_a=0.1,
          D_b=0.1,
          Re_turbulent=100,
          D_Re=0.1,
        zeta2=1) "Same as lossFactors1 but with data for laminar region" annotation (extent=[-100,-82; -80,-62]);
      
      Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_p1(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[40,60; 20,80]);
      Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
        orifice_p1(
        redeclare package Medium = Medium,
        data=lossFactors1) 
        annotation (extent=[-20,60; 0,80]);
      
      annotation (Diagram,
        experiment(StopTime=10, NumberOfIntervals=50000),
        experimentSetupOutput,
        Coordsys(extent=[-100,-200; 200,100]),
        Documentation(info="<html>
<p>
Test whether the different settings of \"from_dp\" and \"use_Re\"
gives the same results for an orifice where the laminar
region is defined and where it is not defined.
</p>
</html>"));
      Modelica_Fluid.Sources.PrescribedBoundary_pTX ambientSource(
                                                  redeclare package Medium = Medium) 
        annotation (extent=[-60,60; -40,80]);
      Modelica.Blocks.Sources.TimeTable p_table(table=[0,0.999e5; 10,1.001e5]) 
        annotation (extent=[-100,60; -80,80]);
      Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_p2(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[40,30; 20,50]);
      Modelica_Fluid.Sources.PrescribedMassFlowRate_TX pump_m1(redeclare 
          package Medium = Medium) 
        annotation (extent=[100,60; 120,80]);
      Modelica.Blocks.Sources.TimeTable m_flow_table(table=[0,-10; 10,10]) 
        annotation (extent=[60,60; 80,80]);
      Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
        orifice_m1(
        redeclare package Medium = Medium,
        data=lossFactors1) 
        annotation (extent=[140,60; 160,80]);
      Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_m1(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[200,60; 180,80]);
      Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
        orifice_p2(
        redeclare package Medium = Medium,
        data=lossFactors1,
        from_dp=false) 
        annotation (extent=[-20,30; 0,50]);
      
      Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_p3(
                                         redeclare package Medium = Medium, p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[40,0; 20,20]);
      Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
        orifice_p3(
        redeclare package Medium = Medium,
        data=lossFactors1,
        use_Re=false) 
        annotation (extent=[-20,0; 0,20]);
      Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_p4(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[40,-30; 20,-10]);
      Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
        orifice_p4(
        redeclare package Medium = Medium,
        data=lossFactors1,
        from_dp=false,
        use_Re=false) 
        annotation (extent=[-20,-30; 0,-10]);
      Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_p5(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[40,-60; 20,-40]);
      Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
        orifice_p5(
        redeclare package Medium = Medium,
        data=lossFactors2) 
        annotation (extent=[-20,-60; 0,-40]);
      Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_p6(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[40,-90; 20,-70]);
      Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
        orifice_p6(
        redeclare package Medium = Medium,
        from_dp=false,
        data=lossFactors2) 
        annotation (extent=[-20,-90; 0,-70]);
      Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_p7(
                                         redeclare package Medium = Medium, p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[40,-120; 20,-100]);
      Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
        orifice_p7(
        redeclare package Medium = Medium,
        use_Re=false,
        data=lossFactors2) 
        annotation (extent=[-20,-120; 0,-100]);
      Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_p8(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[40,-150; 20,-130]);
      Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
        orifice_p8(
        redeclare package Medium = Medium,
        from_dp=false,
        use_Re=false,
        data=lossFactors2) 
        annotation (extent=[-20,-150; 0,-130]);
      Modelica_Fluid.Sources.PrescribedMassFlowRate_TX pump_m2(redeclare 
          package Medium = Medium) 
        annotation (extent=[100,30; 120,50]);
      Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
        orifice_m2(
        redeclare package Medium = Medium,
        data=lossFactors1,
        from_dp=false) 
        annotation (extent=[140,30; 160,50]);
      Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_m2(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[200,30; 180,50]);
      Modelica_Fluid.Sources.PrescribedMassFlowRate_TX pump_m3(redeclare 
          package Medium = Medium) 
        annotation (extent=[100,0; 120,20]);
      Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
        orifice_m3(
        redeclare package Medium = Medium,
        data=lossFactors1,
        use_Re=false) 
        annotation (extent=[140,0; 160,20]);
      Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_m3(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[200,0; 180,20]);
      Modelica_Fluid.Sources.PrescribedMassFlowRate_TX pump_m4(redeclare 
          package Medium = Medium) 
        annotation (extent=[100,-30; 120,-10]);
      Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
        orifice_m4(
        redeclare package Medium = Medium,
        data=lossFactors1,
        from_dp=false,
        use_Re=false) 
        annotation (extent=[140,-30; 160,-10]);
      Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_m4(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[200,-30; 180,-10]);
      Modelica_Fluid.Sources.PrescribedMassFlowRate_TX pump_m5(redeclare 
          package Medium = Medium) 
        annotation (extent=[100,-60; 120,-40]);
      Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
        orifice_m5(
        redeclare package Medium = Medium,
        data=lossFactors2) 
        annotation (extent=[140,-60; 160,-40]);
      Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_m5(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[200,-60; 180,-40]);
      Modelica_Fluid.Sources.PrescribedMassFlowRate_TX pump_m6(redeclare 
          package Medium = Medium) 
        annotation (extent=[100,-90; 120,-70]);
      Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
        orifice_m6(
        redeclare package Medium = Medium,
        from_dp=false,
        data=lossFactors2) 
        annotation (extent=[140,-90; 160,-70]);
      Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_m6(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[200,-90; 180,-70]);
      Modelica_Fluid.Sources.PrescribedMassFlowRate_TX pump_m7(redeclare 
          package Medium = Medium) 
        annotation (extent=[100,-120; 120,-100]);
      Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
        orifice_m7(
        redeclare package Medium = Medium,
        use_Re=false,
        data=lossFactors2) 
        annotation (extent=[140,-120; 160,-100]);
      Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_m7(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[200,-120; 180,-100]);
      Modelica_Fluid.Sources.PrescribedMassFlowRate_TX pump_m8(redeclare 
          package Medium = Medium) 
        annotation (extent=[100,-150; 120,-130]);
      Modelica_Fluid.PressureLosses.BaseClasses.QuadraticTurbulent.BaseModel 
        orifice_m8(
        redeclare package Medium = Medium,
        from_dp=false,
        use_Re=false,
        data=lossFactors2) 
        annotation (extent=[140,-150; 160,-130]);
      Modelica_Fluid.Sources.FixedBoundary_pTX ambientSource_m8(
        redeclare package Medium = Medium,
        p=1.0e5,
        T=Modelica.SIunits.Conversions.from_degC(80)) 
        annotation (extent=[200,-150; 180,-130]);
      inner Modelica_Fluid.Ambient ambient 
                            annotation (extent=[-100,-140; -80,-120]);
    equation 
      connect(orifice_p1.port_b, ambientSource_p1.port) 
        annotation (points=[0,70; 20,70], style(color=69, rgbcolor={0,127,255}));
      connect(ambientSource.port, orifice_p1.port_a) 
                                             annotation (points=[-40,70; -20,70],
          style(color=69, rgbcolor={0,127,255}));
      connect(p_table.y, ambientSource.p_in)  annotation (points=[-79,70; -72,70; -72,76;
            -62,76], style(color=74, rgbcolor={0,0,127}));
      connect(m_flow_table.y, pump_m1.m_flow_in) annotation (points=[81,70; 88,70;
            88,76; 100.7,76], style(color=74, rgbcolor={0,0,127}));
      connect(pump_m1.port, orifice_m1.port_a) annotation (points=[120,70; 140,70],
          style(color=69, rgbcolor={0,127,255}));
      connect(orifice_m1.port_b, ambientSource_m1.port) annotation (points=[160,70; 180,
            70], style(color=69, rgbcolor={0,127,255}));
      connect(ambientSource.port, orifice_p2.port_a) annotation (points=[-40,70; -30,70;
            -30,40; -20,40], style(color=69, rgbcolor={0,127,255}));
      connect(orifice_p2.port_b, ambientSource_p2.port) 
        annotation (points=[0,40; 20,40], style(color=69, rgbcolor={0,127,255}));
      connect(orifice_p3.port_b, ambientSource_p3.port) 
        annotation (points=[0,10; 20,10], style(color=69, rgbcolor={0,127,255}));
      connect(orifice_p4.port_b, ambientSource_p4.port) 
        annotation (points=[0,-20; 20,-20], style(color=69, rgbcolor={0,127,255}));
      connect(ambientSource.port, orifice_p3.port_a) annotation (points=[-40,70; -30,70;
            -30,10; -20,10], style(color=69, rgbcolor={0,127,255}));
      connect(ambientSource.port, orifice_p4.port_a) annotation (points=[-40,70; -30,70;
            -30,-20; -20,-20], style(color=69, rgbcolor={0,127,255}));
      connect(orifice_p5.port_b, ambientSource_p5.port) 
        annotation (points=[0,-50; 20,-50],
                                          style(color=69, rgbcolor={0,127,255}));
      connect(ambientSource.port,orifice_p5. port_a) 
                                             annotation (points=[-40,70; -30,70;
            -30,-50; -20,-50],
          style(color=69, rgbcolor={0,127,255}));
      connect(ambientSource.port, orifice_p6.port_a) annotation (points=[-40,70; -30,70;
            -30,-80; -20,-80], style(color=69, rgbcolor={0,127,255}));
      connect(orifice_p6.port_b, ambientSource_p6.port) 
        annotation (points=[0,-80; 20,-80], style(color=69, rgbcolor={0,127,255}));
      connect(orifice_p7.port_b, ambientSource_p7.port) 
        annotation (points=[0,-110; 20,-110],
                                          style(color=69, rgbcolor={0,127,255}));
      connect(orifice_p8.port_b, ambientSource_p8.port) annotation (points=[0,-140; 20,
            -140], style(color=69, rgbcolor={0,127,255}));
      connect(ambientSource.port, orifice_p7.port_a) annotation (points=[-40,70; -30,70;
            -30,-110; -20,-110], style(color=69, rgbcolor={0,127,255}));
      connect(ambientSource.port, orifice_p8.port_a) annotation (points=[-40,70; -30,70;
            -30,-140; -20,-140], style(color=69, rgbcolor={0,127,255}));
      connect(pump_m2.port, orifice_m2.port_a) annotation (points=[120,40; 140,40],
          style(color=69, rgbcolor={0,127,255}));
      connect(orifice_m2.port_b, ambientSource_m2.port) annotation (points=[160,40; 180,
            40], style(color=69, rgbcolor={0,127,255}));
      connect(pump_m3.port, orifice_m3.port_a) annotation (points=[120,10; 140,10],
          style(color=69, rgbcolor={0,127,255}));
      connect(orifice_m3.port_b, ambientSource_m3.port) annotation (points=[160,10; 180,
            10], style(color=69, rgbcolor={0,127,255}));
      connect(pump_m4.port, orifice_m4.port_a) annotation (points=[120,-20; 140,-20],
          style(color=69, rgbcolor={0,127,255}));
      connect(orifice_m4.port_b, ambientSource_m4.port) annotation (points=[160,-20; 180,
            -20], style(color=69, rgbcolor={0,127,255}));
      connect(m_flow_table.y, pump_m2.m_flow_in) annotation (points=[81,70; 88,70;
            88,46; 100.7,46], style(color=74, rgbcolor={0,0,127}));
      connect(m_flow_table.y, pump_m3.m_flow_in) annotation (points=[81,70; 88,70;
            88,16; 100.7,16], style(color=74, rgbcolor={0,0,127}));
      connect(m_flow_table.y, pump_m4.m_flow_in) annotation (points=[81,70; 88,70;
            88,-14; 100.7,-14], style(color=74, rgbcolor={0,0,127}));
      connect(m_flow_table.y, pump_m5.m_flow_in) annotation (points=[81,70; 88,70;
            88,-44; 100.7,-44], style(color=74, rgbcolor={0,0,127}));
      connect(pump_m5.port, orifice_m5.port_a) annotation (points=[120,-50; 140,-50],
          style(color=69, rgbcolor={0,127,255}));
      connect(orifice_m5.port_b, ambientSource_m5.port) annotation (points=[160,-50; 180,
            -50], style(color=69, rgbcolor={0,127,255}));
      connect(pump_m6.port, orifice_m6.port_a) annotation (points=[120,-80; 140,-80],
          style(color=69, rgbcolor={0,127,255}));
      connect(orifice_m6.port_b, ambientSource_m6.port) annotation (points=[160,-80; 180,
            -80], style(color=69, rgbcolor={0,127,255}));
      connect(pump_m7.port, orifice_m7.port_a) annotation (points=[120,-110; 140,
            -110], style(color=69, rgbcolor={0,127,255}));
      connect(orifice_m7.port_b, ambientSource_m7.port) annotation (points=[160,-110; 180,
            -110], style(color=69, rgbcolor={0,127,255}));
      connect(pump_m8.port, orifice_m8.port_a) annotation (points=[120,-140; 140,
            -140], style(color=69, rgbcolor={0,127,255}));
      connect(orifice_m8.port_b, ambientSource_m8.port) annotation (points=[160,-140; 180,
            -140], style(color=69, rgbcolor={0,127,255}));
      connect(m_flow_table.y, pump_m6.m_flow_in) annotation (points=[81,70; 88,70;
            88,-74; 100.7,-74], style(color=74, rgbcolor={0,0,127}));
      connect(m_flow_table.y, pump_m7.m_flow_in) annotation (points=[81,70; 88,70;
            88,-104; 100.7,-104], style(color=74, rgbcolor={0,0,127}));
      connect(m_flow_table.y, pump_m8.m_flow_in) annotation (points=[81,70; 88,70;
            88,-134; 100.7,-134], style(color=74, rgbcolor={0,0,127}));
    end TestPressureLoss;
    
    model TestPortVolumes 
      extends Modelica.Icons.Example;
      package Medium = Modelica.Media.Water.StandardWater;
      Modelica_Fluid.WorkInProgress.Components.PortVolume PortVolume1(
        V=1e-3,
        use_T_start=false,
        h_start=1e5,
        redeclare package Medium = Medium) 
                     annotation (extent=[-30,-10; -10,10]);
      Modelica_Fluid.Sources.PrescribedMassFlowRate_hX FlowSource1(
        m_flow=1,
        h=2e5,
        redeclare package Medium = Medium) 
                       annotation (extent=[-100,-10; -80,10]);
      Modelica_Fluid.Sources.FixedBoundary_phX Sink1(
                                             p=101325, redeclare package Medium
          =        Medium,
        h=Medium.h_default) 
        annotation (extent=[100,-10; 80,10]);
      Modelica_Fluid.WorkInProgress.Components.PortVolume PortVolume2(
        V=1e-3,
        use_T_start=false,
        h_start=1e5,
        redeclare package Medium = Medium) 
                     annotation (extent=[10,-10; 30,10]);
      Modelica_Fluid.Sensors.TemperatureOnePort Tin(
                                         redeclare package Medium = Medium) 
        annotation (extent=[-60,10; -40,30]);
      Modelica_Fluid.Sensors.TemperatureOnePort Tout(
                                          redeclare package Medium = Medium) 
        annotation (extent=[40,10; 60,30]);
      inner Modelica_Fluid.Ambient ambient 
        annotation (extent=[-100,-100; -80,-80]);
    equation 
      connect(PortVolume1.port, PortVolume2.port) 
        annotation (points=[-20,0; 20,0],style(color=69, rgbcolor={0,127,255}));
      annotation (Diagram);
      connect(PortVolume2.port, Sink1.port) 
        annotation (points=[20,0; 80,0], style(color=69, rgbcolor={0,127,255}));
      connect(PortVolume2.port, Tout.port) annotation (points=[20,0; 50,0; 50,10],
          style(color=69, rgbcolor={0,127,255}));
      connect(FlowSource1.port, PortVolume1.port) 
        annotation (points=[-80,0; -20,0], style(color=69, rgbcolor={0,127,255}));
      connect(FlowSource1.port, Tin.port) annotation (points=[-80,0; -50,0; -50,10],
          style(color=69, rgbcolor={0,127,255}));
    end TestPortVolumes;

    package BatchPlant 
      partial model PartialPump "Base model for centrifugal pumps" 
        import Modelica.SIunits.Conversions.NonSIunits.*;
        import Modelica.Constants.*;
        replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
          "Medium model" annotation(choicesAllMatching=true);
        Medium.BaseProperties fluid(p(start=pin_start),h(start=h_start)) 
          "Fluid properties at the inlet";
        replaceable package SatMedium = 
            Modelica.Media.Interfaces.PartialTwoPhaseMedium 
          "Saturated medium model (required only for NPSH computation)" 
                                                                       annotation(choicesAllMatching=true);
        replaceable function flowCharacteristic = 
            Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics.baseFlow 
          "Head vs. q_flow characteristic at nominal speed and density" 
          annotation(Dialog(group="Characteristics"), choicesAllMatching=true);
        parameter Boolean usePowerCharacteristic = false 
          "Use powerCharacteristic (vs. efficiencyCharacteristic)" 
           annotation(Dialog(group="Characteristics"));
          replaceable function powerCharacteristic = 
            Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics.quadraticPower
            (q_nom={0,0,0},W_nom={0,0,0}) 
          "Power consumption vs. q_flow at nominal speed and density" 
            annotation(Dialog(group="Characteristics", enable = usePowerCharacteristic),
                       choicesAllMatching=true);
        
        replaceable function efficiencyCharacteristic = 
          Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics.constantEfficiency
            ( eta_nom=0.8) extends 
          Modelica_Fluid.Pumps.BaseClasses.PumpCharacteristics.baseEfficiency 
          "Efficiency vs. q_flow at nominal speed and density" 
          annotation(Dialog(group="Characteristics",enable = not usePowerCharacteristic),
                     choicesAllMatching=true);
        parameter AngularVelocity_rpm N_nom = 1500 "Nominal rotational speed" 
          annotation(Dialog(group="Characteristics"));
        parameter Medium.Density d_nom = 1000 "Nominal fluid density" 
          annotation(Dialog(group="Characteristics"));
        parameter Integer Np_nom(min=1) = 1 
          "Nominal number of pumps in parallel";
        parameter Modelica.SIunits.Mass M=0 "Fluid mass inside the pump";
        parameter Boolean checkValve=true "Reverse flow stopped";
        parameter Boolean allowFlowReversal = true 
          "Flow reversal at the ports is allowed by the equations";
        parameter Boolean computeNPSHa=false 
          "Compute NPSH Available at the inlet";
        parameter Medium.AbsolutePressure pin_start 
          "Inlet Pressure Start Value" 
          annotation(Dialog(tab="Initialization"));
        parameter Medium.AbsolutePressure pout_start 
          "Outlet Pressure Start Value" 
          annotation(Dialog(tab="Initialization"));
        parameter Boolean use_T_start = true 
          "Use T_start if true, otherwise h_start" 
          annotation(Dialog(tab = "Initialization"), Evaluate = true);
        parameter Medium.Temperature T_start=
          if use_T_start then 293.15 else Medium.temperature_phX(pin_start,h_start,Medium.reference_X[1:Medium.nXi]) 
          "Start value of temperature" 
          annotation(Dialog(tab = "Initialization", enable = use_T_start));
        parameter Medium.SpecificEnthalpy h_start=
          if use_T_start then Medium.specificEnthalpy_pTX(pin_start, T_start, Medium.reference_X[1:Medium.nXi]) else 1e4 
          "Start value of specific enthalpy" 
          annotation(Dialog(tab = "Initialization", enable = not use_T_start));
        parameter Modelica.SIunits.MassFlowRate m_flow_start=0 
          "Start value of mass flow rate (total)" 
          annotation(Dialog(tab="Initialization"));
        constant Modelica.SIunits.Acceleration g=Modelica.Constants.g_n;
      //  parameter Choices.Init.Options.Temp initOpt=Choices.Init.Options.noInit 
      //    "Initialisation option";
        Modelica_Fluid.Interfaces.FluidPort_a inlet(
          redeclare package Medium = Medium,
          p(start=pin_start),
          m_flow(start=m_flow_start, min=if allowFlowReversal and not checkValve then 
                      -inf else 0)) 
        annotation (extent=[-100,-40; -60,0]);
        Modelica_Fluid.Interfaces.FluidPort_b outlet(
          redeclare package Medium = Medium,
          p(start=pout_start),
          m_flow(start=-m_flow_start, max=if allowFlowReversal and not checkValve then 
                      +inf else 0)) 
        annotation (extent=[40,12; 80,52]);
        Modelica.SIunits.Pressure dp=outlet.p - inlet.p "Pressure increase";
        Modelica.SIunits.Height head=dp/(d*g) "Pump head";
        Medium.Density d "Liquid density at the inlet";
        Medium.SpecificEnthalpy h_out(start=h_start) 
          "Enthalpy of the liquid flowing out of the pump";
        Medium.Temperature Tin "Liquid inlet temperature";
        Modelica.SIunits.MassFlowRate m_flow=inlet.m_flow 
          "Mass flow rate (total)";
        Modelica.SIunits.MassFlowRate m_flow_single=m_flow/Np 
          "Mass flow rate (single pump)";
        Modelica.SIunits.VolumeFlowRate q_flow=m_flow/d 
          "Volume flow rate (total)";
        Modelica.SIunits.VolumeFlowRate q_flow_single=q_flow/Np 
          "Volume flow rate (single pump)";
        AngularVelocity_rpm N "Shaft rotational speed";
        Integer Np(min=1) "Number of pumps in parallel";
        Modelica.SIunits.Power W_single "Power Consumption (single pump)";
        Modelica.SIunits.Power W_tot=W_single*Np "Power Consumption (total)";
        constant Modelica.SIunits.Power W_eps=1e-8 
          "Small coefficient to avoid numerical singularities in efficiency computations";
        Real eta "Global Efficiency";
        Modelica.SIunits.Length NPSHa "Net Positive Suction Head available";
        Medium.AbsolutePressure pv "Saturation pressure of inlet liquid";
        Real s(start = m_flow_start) 
          "Curvilinear abscissa for the flow curve in parametric form (either mass flow rate or head)";
        Modelica.Blocks.Interfaces.IntegerInput in_Np 
          annotation (extent=[16,34; 36,54], rotation=-90);
      protected 
        constant Modelica.SIunits.Height unitHead=1;
        constant Modelica.SIunits.MassFlowRate unitMassFlowRate=1;
      equation 
        // Number of pumps in parallel
        Np = in_Np;
        if cardinality(in_Np)==0 then
          in_Np = Np_nom "Number of pumps selected by parameter";
        end if;
        
        // Flow equations
        if noEvent(s > 0 or (not checkValve)) then
          // Flow characteristics when check valve is open
          // q_flow_single = s;
          q_flow_single = s*unitMassFlowRate/d;
          head = noEvent((((if abs(N) > 1e-6 then N else 1e-6))/N_nom)^2*flowCharacteristic(q_flow_single*N_nom/((if abs(N) > 1e-6 then N else 1e-6))));
        else
          // Flow characteristics when check valve is closed
          head = (N/N_nom)^2*flowCharacteristic(0) - s*unitHead;
          q_flow_single = 0;
        end if;
        
        // Power consumption  
        if usePowerCharacteristic then
          W_single = (N/N_nom)^3*(d/d_nom)*powerCharacteristic(q_flow_single*N_nom/(noEvent(if abs(N) > 1e-6 then N else 1e-6))) 
            "Power consumption (single pump)";
          eta = (dp*q_flow_single)/(W_single + W_eps) "Hydraulic efficiency";
        else
          eta = efficiencyCharacteristic(q_flow_single*N_nom/(noEvent(if abs(N) > 1e-6 then N else 1e-10)));
          W_single = dp*q_flow/eta;
        end if;
        // Fluid properties
        fluid.p = inlet.p;
        fluid.h = inlet.h;
        fluid.Xi = inlet.Xi;
        d = fluid.d;
        Tin = fluid.T;
        
        // Mass and energy balances
        inlet.m_flow + outlet.m_flow = 0 "Mass balance";
        inlet.mXi_flow + outlet.mXi_flow = zeros(Medium.nXi) 
          "Substance mass balance";
        inlet.H_flow=semiLinear(inlet.m_flow,inlet.h,h_out) 
          "Enthalpy flow at the inlet";
        outlet.H_flow=semiLinear(outlet.m_flow,outlet.h,h_out) 
          "Enthalpy flow at the outlet";
        if M > 0 then
          M * der(h_out) = m_flow_single*(inlet.h - outlet.h) + W_single 
            "Dynamic energy balance (density variations neglected)";
        else
          inlet.H_flow + outlet.H_flow + W_single*Np = 0 
            "Static energy balance";
        end if;
        
        // NPSH computations
        if computeNPSHa then
            pv = SatMedium.saturationPressure(fluid.T);
          NPSHa = (inlet.p - pv)/(d*Modelica.Constants.g_n);
        else
          pv = 0;
          NPSHa = 0;
        end if;
      /*
initial equation 
  if initOpt == Choices.Init.Options.noInit then
    // do nothing
  elseif initOpt == Choices.Init.Options.steadyState then
    if ThermalCapacity then
      der(h)=0;
    end if;
  else
    assert(false, "Unsupported initialisation option");
  end if;
*/
        annotation (
          Icon(
            Polygon(points=[-40,-64; -60,-100; 60,-100; 40,-64; -40,-64],
                style(pattern=0, fillColor=74)),
            Ellipse(extent=[-60,40; 60,-80],   style(gradient=3)),
            Polygon(points=[-30,12; -30,-48; 48,-20; -30,12],   style(
                pattern=0,
                gradient=2,
                fillColor=7)),
            Text(extent=[-100,-110; 100,-136], string="%name"),
            Text(extent=[-10,60; 18,40],  string="Np")),
          Diagram,
          Documentation(info="<HTML>
<p>This is the base model for the <tt>Pump</tt> and <tt>
PumpMech</tt> pump models.
<p>The model describes a centrifugal pump, or a group of <tt>Np</tt> identical pumps in parallel. The pump model is based on the theory of kinematic similarity: the pump characteristics are given for nominal operating conditions (rotational speed and fluid density), and then adapted to actual operating condition, according to the similarity equations. 
<p><b>Modelling options</b></p>
<p> The nominal hydraulic characteristic (head vs. volume flow rate) is given by the the replaceable function <tt>flowCharacteristic</tt>. 
<p> The pump energy balance can be specified in two alternative ways:
<ul>
<li><tt>usePowerCharacteristic = false</tt> (default option): the replaceable function <tt>efficiencyCharacteristic</tt> (efficiency vs. volume flow rate in nominal conditions) is used to determine the efficiency, and then the power consumption. The default is a constant efficiency of 0.8.
<li><tt>usePowerCharacteristic = true</tt>: the replaceable function <tt>powerCharacteristic</tt> (power consumption vs. volume flow rate in nominal conditions) is used to determine the power consumption, and then the efficiency.
</ul>
<p>
Several functions are provided in the package <tt>PumpCharacteristics</tt> to specify the characteristics as a function of some operating points at nominal conditions.
<p>Depending on the value of the <tt>checkValve</tt> parameter, the model either supports reverse flow conditions, or includes a built-in check valve to avoid flow reversal.
<p>If the <tt>in_Np</tt> input connector is wired, it provides the number of pumps in parallel; otherwise,  <tt>Np_n</tt> parallel pumps are assumed.</p>
<p>It is possible to take into account the heat capacity of the fluid inside the pump by specifying its mass <tt>M</tt> at nominal conditions; this is necessary to avoid singularities in the computation of the outlet enthalpy in case of zero flow rate. If zero flow rate conditions are always avoided, this dynamic effect can be neglected by leaving the default value <tt>M = 0</tt>, thus avoiding a fast state variable in the model.
<p>If <tt>computeNPSHa = true</tt>, the available net positive suction head is also computed; this requires a two-phase medium model to provide the fluid saturation pressure.
</HTML>",   revisions="<html>
<ul>
<li><i>31 Oct 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
        
      end PartialPump;

      model PumpWithAssertOfCavitation 
        "Centrifugal pump with ideally controlled speed" 
        extends Modelica_Fluid_WorkInProgress.Examples.BatchPlant.PartialPump( final 
            computeNPSHa =                                                                        true);
        import Modelica.SIunits.Conversions.NonSIunits.*;
        parameter AngularVelocity_rpm N_const = N_nom 
          "Constant rotational speed";
        Modelica.Blocks.Interfaces.RealInput N_in "Prescribed rotational speed"
          annotation (extent=[-36,34; -16,54],   rotation=-90);
      equation 
        
         assert(inlet.p >= pv,   " 
    wahrscheinlich ist ein Ventil zu oder ein Tank vor der Pumpe leer.
    ");
        
          N = N_in "Rotational speed";
        if cardinality(N_in)==0 then
          N_in = N_const "Rotational speed provided by parameter";
        end if;
        annotation (
          Icon(
            Text(extent=[-58,58; -30,38], string="n")),
          Diagram,
          Documentation(info="<HTML>
<p>This model describes a centrifugal pump (or a group of <tt>Np</tt> pumps in parallel) with controlled speed, either fixed or provided by an external signal.
<p>The model extends <tt>PartialPump</tt>
<p>If the <tt>N_in</tt> input connector is wired, it provides rotational speed of the pumps (rpm); otherwise, a constant rotational speed equal to <tt>n_const</tt> (which can be different from <tt>N_nom</tt>) is assumed.</p>
</HTML>",   revisions="<html>
<ul>
<li><i>31 Oct 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
      end PumpWithAssertOfCavitation;
    end BatchPlant;
  end Examples;
  
  package FluidStorage 
    
  connector FluidPort_ArrayIcon 
      "Fluid connector with icon suited for an array of FluidPorts" 
      import Modelica_Fluid;
    extends Modelica_Fluid.Interfaces.FluidPort;
    annotation (defaultComponentName="ports",
                Diagram(Rectangle(extent=[-100, 100; 100, -100], style(color=69,
               fillColor=69)), Rectangle(extent=[-100, 100; 100, -100], style(color=16,
               fillColor=69)), Text(extent=[-88, 206; 112, 112], string="%name")),
         Icon(Rectangle(extent=[-100, 100; 100, -100], style(color=69,
              fillColor=69)), Rectangle(extent=[-100, 100; 100, -100], style(color=16,
              fillColor=69))));
  end FluidPort_ArrayIcon;
    
  model TankAttachment "Equations to attach pipe at tank" 
        replaceable package Medium = 
        Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching=true);
      
      Modelica_Fluid.Interfaces.FluidPort_a port(redeclare package Medium = Medium) 
      annotation (extent=[-10,-112; 10,-92],    rotation=90);
     // Real mXi_flow;
      parameter Boolean onlyInFlow = false 
        "= true, if flow only into the tank (e.g. top ports)" 
                                                            annotation(Evaluate=true);
      parameter Modelica.SIunits.Diameter pipeDiameter=0.0 
        "Inner (hydraulic) pipe diameter"                                  annotation(Dialog(enable=not onlyInFlow));
      parameter Modelica.SIunits.Height pipeHeight "Height of pipe";
      parameter Medium.SpecificEnthalpy h_start 
        "Start value of specific enthalpy (used as enthalpy for topPorts if back flow)";
      parameter Medium.MassFraction X_start[Medium.nX] 
        "Start value of mass fractions m_i/m";
      parameter Modelica.SIunits.Height level_start "Start value of tank level";
      
      parameter Modelica.SIunits.AbsolutePressure p_ambient 
        "Tank surface pressure"                                   annotation(Dialog);
      input Modelica.SIunits.Height level "Actual tank level" 
                                                annotation(Dialog);
      input Medium.SpecificEnthalpy h 
        "Actual specific enthalpy of fluid in tank"               annotation(Dialog);
      input Medium.Density d "Actual specific density of fluid in tank" 
                                                        annotation(Dialog);
      input Medium.MassFraction Xi[Medium.nXi] 
        "Actual mass fractions of fluid in tank"                  annotation(Dialog);
      parameter Real k_small(min=0) = 1e-5 
        "Small regularization range if tank level is below bottom_height or side_height; k_small = 0 gives ideal switch"
                annotation(Evaluate=true);
      parameter Real s_start = 0;
      
      output Medium.EnthalpyFlowRate H_flow 
        "= port.H_flow (used to transform vector of connectors in vector of Real numbers)";
      output Medium.MassFlowRate m_flow 
        "= port.m_flow (used to transform vector of connectors in vector of Real numbers)";
      output Medium.MassFlowRate mXi_flow[Medium.nXi] 
        "= port.mXi_flow (used to transform vector of connectors in vector of Real numbers)";
      
    annotation (Documentation(info="<html>
<p>
This component contains the equations that attach the pipe
to the tank. The main reason to introduce this component is
that Dymola has currently limitations for connector arrays
when the dimension is zero. Without this utility component
it would not be possible to set, e.g., n_topPorts to zero.
</p>
 
 
</html>"), Icon(Rectangle(extent=[-100,0; 100,-100], style(
            color=3,
            rgbcolor={0,0,255},
            fillColor=7,
            rgbfillColor={255,255,255})), Text(
          extent=[-122,48; 132,6],
          style(
            color=3,
            rgbcolor={0,0,255},
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1),
          string="%name")));
    Modelica_Fluid.Ambient ambient;
    protected 
    parameter Modelica.SIunits.Area pipeArea=Modelica.Constants.pi*(
          pipeDiameter/2)^2;
    parameter Medium.MassFlowRate m_flow_nominal = 1 
        "Nominal mass flow rate used for scaling (has only an effect if k_small > 0)";
    parameter Medium.AbsolutePressure p_nominal = p_ambient;
    Modelica.SIunits.Length aboveLevel=level - pipeHeight;
    Boolean m_flow_out(start=true,fixed=true) "true= massflow out of tank";
    Real s(start=s_start) 
        "path parameter of parameterized curve description (either m_flow/m_flow_nominal or (port.p-p_ambient)/p_ambient)";
  equation 
    H_flow = port.H_flow;
    m_flow = port.m_flow;
    mXi_flow = port.mXi_flow;
      
    if onlyInFlow then
       m_flow_out = false "Dummy value in this if branch";
       port.p = p_ambient;
       /* flow should never out of the port. However, if this occurs in a 
        small time interval (e.g. during initialization), the start values of
        h and X are provided, since otherwise there is a singular
        system 
     */
       port.H_flow = semiLinear(port.m_flow, port.h, h_start);
       port.mXi_flow = semiLinear(port.m_flow, port.Xi, X_start[1:Medium.nXi]);
       assert(port.m_flow > -1e-6, "Mass flows out of tank via topPort. This indicates a wrong model");
       s = 0;
    else
       port.H_flow = semiLinear(port.m_flow, port.h, h);
       port.mXi_flow = semiLinear(port.m_flow, port.Xi, Xi);
        
  /* Original equations from Poschlad/Remelhe:
*/
       s = 0;
       m_flow_out = (pre(m_flow_out) and not port.p>p_ambient) or (port.m_flow < -1e-6);
        
       if (aboveLevel > 0) then
         port.p = aboveLevel*ambient.g*d + p_ambient - smooth(2,noEvent(if m_flow < 0 then m_flow^2/(2*d*pipeArea^2) else 0));
       else
         if pre(m_flow_out) then
            m_flow = 0;
         else
            port.p = p_ambient;
         end if;
       end if;
        
  /* Martin Otter: The following equations are a declarative form 
   (parameterized curve description) of the above equations and
   should theoretically work better. However, some examples with
   IF97 water fail, whereas the above works. Therefore, not used. 
       Add the following text to OpenTank, once the initialization
   with this solution works for Modelica_Fluid.Examples.Tanks.ThreeOpenTanks:
 
OpenTank:
<p>
The situation when the tank level is below bottom_heights[i] or side_heights[i]
is handeled properly. Details are described
<a href="Modelica:Modelica_Fluid.Utilities.TankAttachment">here</a>
</p> 
 
TankAttachment:
<p>
If a bottom or side connector is above the actual tank level, the
following characteristic is used to compute the mass flow rate port.m_flow
from the connector to the tank and the absolute pressure port.p
in the port:
</p>
 
<img src="../Images/Components/Tank_PipeAboveTankLevel.png">   
 
 
     m_flow_out = s <= 0;
     if aboveLevel >= 0 then
        m_flow = m_flow_nominal*s "equation to compute s, which is a dummy in this branch";
        port.p - p_ambient = aboveLevel*fluidOptions.g*d  -
                             smooth(2,if m_flow_out then s*abs(s)*m_flow_nominal^2/(2*d*pipeArea^2) else k_small*m_flow_nominal*s);
     else
        m_flow = (if m_flow_out then k_small*p_nominal else m_flow_nominal)*s;
        port.p - p_ambient = (if m_flow_out then p_nominal else k_small*m_flow_nominal)*s;
     end if;
*/
        
    end if;
      
    /*
  More precise equations (introduce them later; need to transform
  them from energy balance form 1 to form 2):
 
  Momentum balance:
  (integrated momentum equation for frictionless fluid with density that is
   independent of the level, i.e., the unsteady Bernoulli equation for incompressible fluid)
  v_level = der(level);
  v = -port.m_flow/(rho*A_outlet);
  level*der(v_level) + (v^2 - v_level^2)/2 - g*level + (p - p_ambient)/rho = 0; or
  rho*level*der(v_level) + rho*(v^2 - v_level^2)/2 - rho*g*level + (p - p_ambient) = 0;
 
  Energy balance:
  Potential energy: E_pot = integ(dm*g*s)
                          = g*integ(rho*A*s*ds)
                          = g*rho*A*z^2/2
  Kinetic energy  : E_kin = integ(dm*v^2/2)
                          = integ(rho*A*v^2/2*ds)
                          = rho*A*v^2/2*integ(ds)
                          = rho*A*v^2/2*z
                          = M*v^2/2
  E = U + M*g*z/2 + M*v_level^2/2
  der(E) = port.H_flow + port.m_flow*v^2/2 - p_ambient*area*der(level)
*/
      
  end TankAttachment;
  end FluidStorage;
  
  package Interfaces 
    partial model PartialTwoPortTransport 
      "Partial element transporting fluid between two ports without storing mass or energy" 
      import Modelica.SIunits.*;
      import Modelica.Constants.*;
      import Modelica_Fluid;
      replaceable package Medium = 
        Modelica.Media.Interfaces.PartialMedium "Medium in the component"  annotation (
          choicesAllMatching =                                                                            true);
      parameter Boolean allowFlowReversal = true 
        "Flow reversal at the ports is allowed by the equations"  annotation(Dialog(tab="Advanced"));
      parameter Modelica_Fluid.Types.FlowDirection.Temp flowDirection=
                       Modelica_Fluid.Types.FlowDirection.Bidirectional 
        "Unidirectional (port_a -> port_b) or bidirectional flow component" 
         annotation(Dialog(tab="Advanced"));
      
      Modelica_Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = 
            Medium, m_flow(min=if allowFlowReversal then -inf else 0)) 
        annotation (extent=[-120, -10; -100, 10]);
      Modelica_Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = 
            Medium, m_flow(max=if allowFlowReversal then +inf else 0)) 
        annotation (extent=[120, -10; 100, 10]);
      Medium.BaseProperties medium_a "Medium properties in port_a";
      Medium.BaseProperties medium_b "Medium properties in port_b";
      annotation (
        Coordsys(grid=[1, 1], component=[20, 20]),
        Diagram,
        Documentation(info="<html>
<p>
This component transports fluid between its two ports, without
storing mass or energy. Reversal and zero mass flow rate is taken
care of, for details see definition of built-in operator semiLinear().
<p>
When using this partial component, an equation for the momentum
balance has to be added by specifying a relationship
between the pressure drop <tt>dp</tt> and the mass flow rate <tt>m_flow</tt>.
</p>
</html>"));
    equation 
      // Properties in the ports
      port_a.p   = medium_a.p;
      port_a.h   = medium_a.h;
      port_a.Xi = medium_a.Xi;
      port_b.p   = medium_b.p;
      port_b.h   = medium_b.h;
      port_b.Xi = medium_b.Xi;
      
      /* Handle reverse and zero flow */
      port_a.H_flow   = semiLinear(port_a.m_flow, port_a.h,  port_b.h);
      port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.Xi, port_b.Xi);
      
      /* Energy, mass and substance mass balance */
      port_a.H_flow + port_b.H_flow = 0;
      port_a.m_flow + port_b.m_flow = 0;
      port_a.mXi_flow + port_b.mXi_flow = zeros(Medium.nXi);
    end PartialTwoPortTransport;
    
    partial model PartialPressureLoss 
      "Generic pressure loss with constant turbulent loss factors" 
      import Modelica_Fluid.Utilities.regRoot2;
      import Modelica_Fluid.Utilities.regSquare2;
      
      parameter Modelica_Fluid_WorkInProgress.Utilities.PressureLossFactors 
        lossFactors "Loss factors for both flow directions";
      parameter Boolean from_dp=true 
        " = true, use m_flow = f(dp) else dp = f(m_flow)" 
        annotation (Evaluate=true, Dialog(tab="Advanced"));
      parameter Boolean use_Re = true 
        "= true, if turbulent region is defined by Re, otherwise by p_small or m_flow_small"
        annotation(Evaluate=true, Dialog(tab="Advanced"));
      parameter Modelica.SIunits.AbsolutePressure dp_small=1 
        "Turbulent flow if |dp| >= dp_small" 
        annotation(Dialog(tab="Advanced", enable=not use_Re and from_dp));
      parameter Modelica.SIunits.MassFlowRate m_flow_small=0.001 
        "Turbulent flow if |m_flow| >= m_flow_small" 
        annotation(Dialog(tab="Advanced", enable=not use_Re and not from_dp));
      Modelica.SIunits.MassFlowRate m_flow(start=0) 
        "Mass flow rate from port_a to port_b";
      Modelica.SIunits.Pressure dp(start=0) 
        "Pressure loss due to pressure loss component";
      input Modelica.SIunits.Density d_a "Density at port_a" 
                                                   annotation(Hide=true);
      input Modelica.SIunits.Density d_b "Density at port_b" 
                                                   annotation(Hide=true);
      input Modelica.SIunits.DynamicViscosity eta_a 
        "Dynamic viscosity at port_a, if use_Re = true (otherwise not used)";
      input Modelica.SIunits.DynamicViscosity eta_b 
        "Dynamic viscosity at port_b, if use_Re = true (otherwise not used)";
      final Modelica.SIunits.ReynoldsNumber Re=noEvent(abs(m_flow))*(4/pi)/(
          lossFactors.D_Re*eta) if use_Re 
        "Reynolds number at smallest pipe diameter";
      Modelica.SIunits.AbsolutePressure dp_turbulent 
        "The turbulent region is: |dp| >= dp_turbulent";
      Modelica.SIunits.MassFlowRate m_flow_turbulent 
        "The turbulent region is: |m_flow| >= m_flow_turbulent";
      annotation (
        Diagram,
        Icon,
        Documentation(info="<html>
</html>"));
      
    protected 
      constant Real pi=Modelica.Constants.pi annotation(Hide=true);
      parameter Modelica.SIunits.Diameter LD_a=lossFactors.D_a 
                                                  annotation(Hide=true);
      parameter Modelica.SIunits.Diameter LD_b=lossFactors.D_b 
                                                  annotation(Hide=true);
      parameter Real k0=2*lossFactors.c0/(pi*lossFactors.D_Re^3);
      parameter Real k1=if lossFactors.zeta1_at_a then 8*lossFactors.zeta1/(pi*LD_a^2)^2 else 
                                                       8*lossFactors.zeta1/(pi*LD_b^2)^2;
      parameter Real k2=if lossFactors.zeta2_at_a then 8*lossFactors.zeta2/(pi*LD_a^2)^2 else 
                                                       8*lossFactors.zeta2/(pi*LD_b^2)^2;
      parameter Boolean withLaminar = lossFactors.zetaLaminarKnown annotation(Evaluate=true,Hide=true);
      Real yd0 
        "Derivative of dp=dp(m_flow) or m_flow=m_flow(dp) at zero, if use_Re and lossFactors.zetaLaminarKnown";
      Modelica.SIunits.DynamicViscosity eta=if use_Re then (eta_a + eta_b)/2 else 
                0;
      
    equation 
    /*
Turbulent region:
   Re = m_flow*(4/pi)/(D_Re*eta)  
   dp = 0.5*zeta*d*v*|v|
      = 0.5*zeta*d*1/(d*A)^2 * m_flow * |m_flow|
      = 0.5*zeta/A^2 *1/d * m_flow * |m_flow|
      = k/d * m_flow * |m_flow| 
   k  = 0.5*zeta/A^2
      = 0.5*zeta/(pi*(D/2)^2)^2
      = 8*zeta/(pi*D^2)^2
   m_flow_turbulent = (pi/4)*D_Re*eta*Re_turbulent
   dp_turbulent     =  k/d *(D_Re*eta*pi/4)^2 * Re_turbulent^2
_
   The start of the turbulent region is computed with mean values
   of dynamic viscosity eta and density rho. Otherwise, one has
   to introduce different "delta" values for both flow directions.
   In order to simplify the approach, only one delta is used.  
_
Laminar region:
   dp = 0.5*zeta/(A^2*d) * m_flow * |m_flow|
      = 0.5 * c0/(|m_flow|*(4/pi)/(D_Re*eta)) / ((pi*(D_Re/2)^2)^2*d) * m_flow*|m_flow|
      = 0.5 * c0*(pi/4)*(D_Re*eta) * 16/(pi^2*D_Re^4*d) * m_flow*|m_flow|
      = 2*c0/(pi*D_Re^3) * eta/d * m_flow
      = k0 * eta/d * m_flow
   k0 = 2*c0/(pi*D_Re^3)
_
   In order that the derivative of dp=f(m_flow) is continuous 
   at m_flow=0, the mean values of eta and d are used in the
   laminar region: eta/d = (eta_a + eta_b)/(d_a + d_b)
   If lossFactors.zetaLaminarKnown = false then eta_a and eta_b are potentially zero
   (because dummy values) and therefore the division is only performed
   if zetaLaminarKnown = true.
*/
       if from_dp then
          dp_turbulent = if use_Re then (k1+k2)/(d_a+d_b)*(eta*lossFactors.D_Re*pi/4)^2
                                       *lossFactors.Re_turbulent^2 else dp_small;
          m_flow_turbulent = regRoot2(dp_turbulent, dp_turbulent, d_a/k1, d_b/k2) 
          "for information purposes";
          yd0 = if use_Re and withLaminar then (d_a+d_b)/(k0*(eta_a+eta_b)) else 0;
          m_flow = regRoot2(dp, dp_turbulent, d_a/k1, d_b/k2, use_Re and withLaminar, yd0);
       else
          m_flow_turbulent = if use_Re then (pi/4)*lossFactors.D_Re*eta*lossFactors.Re_turbulent else m_flow_small;
          dp_turbulent = regSquare2(m_flow_turbulent, m_flow_turbulent, k1/d_a, k2/d_b) 
          "for information purposes";
          yd0 = if use_Re and withLaminar then k0*(eta_a+eta_b)/(d_a+d_b) else 0;
          dp = regSquare2(m_flow, m_flow_turbulent, k1/d_a, k2/d_b, use_Re and withLaminar, yd0);
       end if;
      
    end PartialPressureLoss;
    
    model PressureLossWithoutIcon 
      "Generic pressure loss component with constant turbulent loss factors and without an icon" 
      extends Modelica_Fluid_WorkInProgress.Interfaces.PartialTwoPortTransport;
      extends Modelica_Fluid_WorkInProgress.Interfaces.PartialPressureLoss(
        m_flow=port_a.m_flow,
        dp=port_a.p - port_b.p,
        d_a=medium_a.d,
        d_b=medium_b.d,
        eta_a=if use_Re then Medium.dynamicViscosity(medium_a.state) else 0,
        eta_b=if use_Re then Medium.dynamicViscosity(medium_b.state) else 0);
      annotation (
        Diagram,
        Icon,
        Documentation(info="<html>
<p>
This model computes the pressure loss of a pipe
segment (orifice, bending etc.) with a minimum amount of data
provided via parameter <b>lossFactors</b>.
If available, data should be provided for <b>both flow directions</b>,
i.e., flow from port_a to port_b and from port_b to port_a, 
as well as for the <b>laminar</b> and the <b>turbulent</b> region.
It is also an option to provide the loss factor <b>only</b> for the
<b>turbulent</b> region for a flow from port_a to port_b.
</p>
<p>
The following equations are used:
</p>
<pre>   &Delta;p = 0.5*&zeta;*&rho;*v*|v|
      = 0.5*&zeta;/A^2 * (1/&rho;) * m_flow*|m_flow|
        Re = |v|*D*&rho;/&eta;
</pre>
<table border=1 cellspacing=0 cellpadding=2>
<tr><td><b>flow type</b></td>
    <td><b>&zeta;</b> = </td>
    <td><b>flow region</b></td></tr>
<tr><td>turbulent</td>
    <td><b>zeta1</b> = const.</td>
    <td>Re &ge;  Re_turbulent, v &ge; 0</td></tr>
<tr><td></td>
    <td><b>zeta2</b> = const.</td>
    <td>Re &ge; Re_turbulent, v &lt; 0</td></tr>
<tr><td>laminar</td>
    <td><b>c0</b>/Re</td>
    <td>both flow directions, Re small; c0 = const.</td></tr>
</table>
<p>
where
</p>
<ul>
<li> &Delta;p is the pressure drop: &Delta;p = port_a.p - port_b.p</li>
<li> v is the mean velocity.</li>
<li> &rho; is the density.</li>
<li> &zeta; is the loss factor that depends on the geometry of
     the pipe. In the turbulent flow regime, it is assumed that
     &zeta; is constant and is given by \"zeta1\" and
     \"zeta2\" depending on the flow direction.<br>
     When the Reynolds number Re is below \"Re_turbulent\", the
     flow is laminar for small flow velocities. For higher 
     velocities there is a transition region from 
     laminar to turbulent flow. The loss factor for
     laminar flow at small velocities is defined by the often occuring
     approximation c0/Re. If c0 is different for the two
     flow directions, the mean value has to be used 
     (c0 = (c0_ab + c0_ba)/2).<li>
<li> The equation \"&Delta;p = 0.5*&zeta;*&rho;*v*|v|\" is either with
     respect to port_a or to port_b, depending on the definition
     of the particular loss factor &zeta; (in some references loss
     factors are defined with respect to port_a, in other references
     with respect to port_b).</li>
 
<li> Re = |v|*D_Re*&rho;/&eta; = |m_flow|*D_Re/(A_Re*&eta;) 
     is the Reynolds number at the smallest cross
     section area. This is often at port_a or at port_b, but can
     also be between the two ports. In the record, the diameter
     D_Re of this smallest cross section area has to be provided, as
     well, as Re_turbulent, the absolute value of the 
     Reynolds number at which
     the turbulent flow starts. If Re_turbulent is different for
     the two flow directions, use the smaller value as Re_turbulent.</li>
<li> D is the diameter of the pipe. If the pipe has not a 
     circular cross section, D = 4*A/P, where A is the cross section
     area and P is the wetted perimeter.</li>
<li> A is the cross section area with A = &pi;(D/2)^2.
<li> &eta; is the dynamic viscosity.</li>
</ul>
<p>
The laminar and the transition region is usually of
not much technical interest because the operating point is
mostly in the turbulent regime. For simplification and for
numercial reasons, this whole region is described by two
polynomials of third order, one polynomial for m_flow &ge; 0 
and one for m_flow &lt; 0. The polynomials start at 
Re = |m_flow|*4/(&pi;*D_Re*&eta;), where D_Re is the
smallest diameter between port_a and port_b.
The common derivative
of the two polynomials at Re = 0 is
computed from the equation \"c0/Re\". Note, the pressure drop
equation above in the laminar region is always defined
with respect to the smallest diameter D_Re.
</p>
<p>
If no data for c0 is available, the derivative at Re = 0 is computed in such
a way, that the second derivatives of the two polynomials
are identical at Re = 0. The polynomials are constructed, such that
they smoothly touch the characteristic curves in the turbulent
regions. The whole characteristic is therefore <b>continuous</b>
and has a <b>finite</b>, <b>continuous first derivative everywhere</b>.
In some cases, the constructed polynomials would \"vibrate\". This is 
avoided by reducing the derivative at Re=0 in such a way that
the polynomials are guaranteed to be monotonically increasing.
The used sufficient criteria for monotonicity follows from:
</p>
 
<dl>
<dt> Fritsch F.N. and Carlson R.E. (1980):</dt>
<dd> <b>Monotone piecewise cubic interpolation</b>.
     SIAM J. Numerc. Anal., Vol. 17, No. 2, April 1980, pp. 238-246</dd>
</dl>
</html>"));
    end PressureLossWithoutIcon;
    
    partial model PartialTwoPortTransportWithDz 
      "Partial element transporting fluid between two ports without storing mass or energy" 
      import Modelica.SIunits.*;
      import Modelica.Constants.*;
      replaceable package Medium = 
          Modelica.Media.Interfaces.PartialTwoPhaseMedium                          extends 
        Modelica.Media.Interfaces.PartialTwoPhaseMedium 
        "Medium in the component"                                          annotation (
          choicesAllMatching =                                                                            true);
      parameter Boolean allowFlowReversal = true 
        "Flow reversal at the ports is allowed by the equations";
      
      Modelica_Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium = 
            Medium, m_flow(min=if allowFlowReversal then -inf else 0)) 
        annotation (extent=[-120, -10; -100, 10]);
      Modelica_Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium = 
            Medium, m_flow(max=if allowFlowReversal then +inf else 0)) 
        annotation (extent=[120, -10; 100, 10]);
      Medium.BaseProperties medium_a "Medium properties in port_a";
      Medium.BaseProperties medium_b "Medium properties in port_b";
      Medium.MassFlowRate m_flow 
        "Mass flow rate from port_a to port_b (m_flow > 0 is design flow direction)";
      Pressure dp(start=0) "Pressure difference between port_a and port_b";
      Real dz_in=0 "dz=hb-ha difference of height";
      constant Modelica.SIunits.Acceleration g=Modelica.Constants.g_n;
      parameter Medium.AbsolutePressure p_ambient=101325 "ambient pressure";
      Medium.SaturationProperties sat 
        "State vector to compute saturation properties";
      Medium.Density rho_l=Medium.bubbleDensity(sat) "density in liquid phase";
      Medium.Density rho_v=Medium.dewDensity(sat) "density in liquid phase";
      Integer help;
      Boolean liquid( start=true);
      annotation (
        Coordsys(grid=[1, 1], component=[20, 20]),
        Diagram,
        Documentation(info="<html>
<p>
This component transports fluid between its two ports, without
storing mass or energy. Reversal and zero mass flow rate is taken
care of, for details see definition of built-in operator semiLinear().
<p>
When using this partial component, an equation for the momentum
balance has to be added by specifying a relationship
between the pressure drop <tt>dp</tt> and the mass flow rate <tt>m_flow</tt>.
</p>
</html>"));
    equation 
      // Properties in the ports
      port_a.p   = medium_a.p;
      port_a.h   = medium_a.h;
      port_a.Xi = medium_a.Xi;
      port_b.p   = medium_b.p;
      port_b.h   = medium_b.h;
      port_b.Xi = medium_b.Xi;
      
    if m_flow > 0 then
      if pre(liquid) then
        if medium_a.d < 1*rho_v+0*rho_l then
          liquid = false;
          help = 1;
        else
          liquid = true;
          help = 2;
        end if;
      else
        if medium_a.d > 1*rho_l+0*rho_v then
          liquid = true;
          help = 3;
        else
          liquid = false;
          help = 4;
        end if;
      end if;
    else
      if pre(liquid) then
        if medium_b.d < 1*rho_v+0*rho_l then
          liquid = false;
          help = 5;
        else
          liquid = true;
          help = 6;
        end if;
      else
        if medium_b.d > 1*rho_l+0*rho_v then
          liquid = true;
          help = 7;
        else
          liquid = false;
          help = 8;
        end if;
      end if;
    end if;
      
      sat.psat = p_ambient;
      sat.Tsat = Medium.saturationTemperature(p_ambient);
      
      /* Handle reverse and zero flow */
      port_a.H_flow   = semiLinear(port_a.m_flow, port_a.h,  port_b.h);
      port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.Xi, port_b.Xi);
      
      /* Energy, mass and substance mass balance */
      port_a.H_flow + port_b.H_flow = 0;
      port_a.m_flow + port_b.m_flow = 0;
      port_a.mXi_flow + port_b.mXi_flow = zeros(Medium.nXi);
      
      // Design direction of mass flow rate
      m_flow = port_a.m_flow;
      
      // Pressure difference between ports
      
      dp = port_a.p - port_b.p - dz_in*g*(if pre(liquid) then rho_l else rho_v);
      
    end PartialTwoPortTransportWithDz;
    
    partial model PartialInitializationParameters 
      "Define parameter menu to initialize medium in component that has one medium model" 
      import Modelica_Fluid.Types;
      
      replaceable package Medium = 
        Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
          annotation (choicesAllMatching = true);
      parameter Types.Init.Temp initType=
                Types.Init.NoInit "Initialization option" 
        annotation(Evaluate=true, Dialog(tab = "Initialization"));
      parameter Medium.AbsolutePressure p_start = Medium.p_default 
        "Start value of pressure" 
        annotation(Dialog(tab = "Initialization"));
      parameter Boolean use_T_start = true 
        "= true, use T_start, otherwise h_start" 
        annotation(Dialog(tab = "Initialization"), Evaluate=true);
      parameter Medium.Temperature T_start=
        if use_T_start then Medium.T_default else Medium.temperature_phX(p_start,h_start,X_start) 
        "Start value of temperature" 
        annotation(Dialog(tab = "Initialization", enable = use_T_start));
      parameter Medium.SpecificEnthalpy h_start=
        if use_T_start then Medium.specificEnthalpy_pTX(p_start, T_start, X_start) else Medium.h_default 
        "Start value of specific enthalpy" 
        annotation(Dialog(tab = "Initialization", enable = not use_T_start));
      parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
        "Start value of mass fractions m_i/m" 
        annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
      annotation (uses(Modelica(version="2.2.2"), Modelica_Fluid(version=
                "1.0 Beta 2")));
    end PartialInitializationParameters;
  end Interfaces;
  
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
    
    package OrificeCharacteristics "Functions for valve characteristics" 
      
      annotation (preferedView="info",
    Documentation(info="<html>
<p>
This package provides functions that compute a 
<b>constant pressure loss factor</b> for 
<b>turbulent flow</b> based solely on geometrical
information of the corresponding component.
All functions have to be derived from function 
<b>partialTurbulentLossFactor</b>.
</p>
<p>
The turbulent pressure loss factor is defined
according to the following equations:
</p>
<pre>
  dp       = port_a.p - port_b.p     \"pressure drop from port_a to port_b\"
  m_flow_a = port_a.m_flow           \"mass flow rate into port_a\"
  d_a      = d_a(port_a.p, port_a.h) \"density at port_a\" 
  v_a      = m_flow_a/(d_a*A_a)      \"velocity at port_a\"
&nbsp;
  // It is assumed that m_flow_a > 0
  dp = 0.5*d_a*zeta_a*v_a^2
     = 0.5*d_a*zeta_a*m_flow_a^2/(d_a*A_a)^2
     = (0.5*zeta_a/A_a^2) * (1/d_a) * m_flow_a^2
</pre>
<p>
When the flow is reversed, the pressure loss equation is defined
in the other way round, using the quantities at port_b.
This gives the following overall description:
</p>
<pre>
  dp = <b>if</b> m_flow_a &gt; 0 <b>then</b>
           k_a * 1/d_a * m_flow_a^2
       <b>else</b>
          -k_b * 1/d_b * m_flow_b^2
&nbsp;
       k_a = 0.5*zeta_a/A_a^2
       k_b = 0.5*zeta_b/A_b^2
</pre>
<p>
This functions returns the pressure loss factors k_a and k_b
as a vector with two elements:
</p>
<pre>
  k[1]: loss factor k_a if m_flow_a &gt; 0
  k[2]: loss factor k_b if m_flow_a &lt; 0
</pre>
<p>
If a pressure loss factor is only known for one flow
direction, m_flow_a &gt; 0, the component can only be
used, if the flow is most of the time in this direction.
One might define k[2]=k[1], if only for short time periods
also a reversal flow appears and then neglect the error
for the flow reversal and at small mass flow rates
when the flow is laminar.
</p>
</html>"));
      
      model turbulentLossFactor "Constant loss factors for turbulent flow" 
        parameter Real k_a "Loss factor if fluid flows from port_a to port_b";
        parameter Real k_b=k_a 
          "Loss factor if fluid flows from port_b to port_a";
      end turbulentLossFactor;
      
      model suddenExpansion "Suddenly expanding area" 
            parameter Modelica.SIunits.Area A_a "Area at port_a";
        parameter Modelica.SIunits.Area A_b 
          "Area at port_b (A_b > A_a required)";
        extends turbulentLossFactor(
                final k_a = (1 - A_a/A_b)^2/A_a^2,
                final k_b = 0.5*(1 - A_a/A_b)^0.75/A_a^2);
        annotation (Icon(
            Rectangle(extent=[-100,40; 0,-40], style(
                color=7,
                rgbcolor={255,255,255},
                fillColor=7,
                rgbfillColor={255,255,255})),
            Rectangle(extent=[0,100; 100,-100], style(
                color=7,
                rgbcolor={255,255,255},
                fillColor=7,
                rgbfillColor={255,255,255})),
            Line(points=[0,40; -100,40; -100,-40; 0,-40; 0,-100; 100,-100; 100,100; 0,
                  100; 0,40], style(
                color=0,
                rgbcolor={0,0,0},
                fillColor=7,
                rgbfillColor={255,255,255},
                fillPattern=1)),
            Line(points=[-76,4; 34,4], style(
                color=3,
                rgbcolor={0,0,255},
                fillColor=7,
                rgbfillColor={255,255,255},
                fillPattern=1)),
            Polygon(points=[34,16; 34,-10; 74,4; 34,16], style(
                color=3,
                rgbcolor={0,0,255},
                fillColor=3,
                rgbfillColor={0,0,255},
                fillPattern=1))), Diagram(
            Line(points=[0,40; -100,40; -100,-40; 0,-40; 0,-100; 100,-100; 100,100; 0,
                  100; 0,40], style(
                color=0,
                rgbcolor={0,0,0},
                fillColor=7,
                rgbfillColor={255,255,255},
                fillPattern=1)),
            Rectangle(extent=[-100,40; 0,-40], style(
                color=7,
                rgbcolor={255,255,255},
                fillColor=7,
                rgbfillColor={255,255,255})),
            Rectangle(extent=[0,100; 100,-100], style(
                color=7,
                rgbcolor={255,255,255},
                fillColor=7,
                rgbfillColor={255,255,255})),
            Line(points=[0,40; -100,40; -100,-40; 0,-40; 0,-100; 100,-100; 100,100; 0,
                  100; 0,40], style(
                color=0,
                rgbcolor={0,0,0},
                fillColor=7,
                rgbfillColor={255,255,255},
                fillPattern=1)),
            Line(points=[-60,-40; -60,40], style(
                color=3,
                rgbcolor={0,0,255},
                arrow=3,
                fillColor=3,
                rgbfillColor={0,0,255},
                fillPattern=1)),
            Text(
              extent=[-50,16; -26,-10],
              style(
                color=3,
                rgbcolor={0,0,255},
                arrow=3,
                fillColor=3,
                rgbfillColor={0,0,255},
                fillPattern=1),
              string="A_a"),
            Line(points=[34,-100; 34,100], style(
                color=3,
                rgbcolor={0,0,255},
                arrow=3,
                fillColor=3,
                rgbfillColor={0,0,255},
                fillPattern=1)),
            Text(
              extent=[54,16; 78,-10],
              style(
                color=3,
                rgbcolor={0,0,255},
                arrow=3,
                fillColor=3,
                rgbfillColor={0,0,255},
                fillPattern=1),
              string="A_b")));
      equation 
        
      end suddenExpansion;
      
    end OrificeCharacteristics;
  end Types;
  
  package Utilities 
  model TankAttachment "Equations to attach pipe at tank" 
         replaceable package Medium = Modelica.Media.Interfaces.PartialMedium 
        "Medium in the component" annotation (choicesAllMatching=true);
      
      Modelica_Fluid.Interfaces.FluidPort_a port(redeclare package Medium = Medium) 
      annotation (extent=[-10,-112; 10,-92],    rotation=90);
     // Real mXi_flow;
      parameter Boolean onlyInFlow = false 
        "= true, if flow only into the tank (e.g. top ports)" 
                                                            annotation(Evaluate=true);
      parameter Modelica.SIunits.Diameter pipeDiameter=0.0 
        "Inner (hydraulic) pipe diameter"                                  annotation(Dialog(enable=not onlyInFlow));
      parameter Modelica.SIunits.Height pipeHeight "Height of pipe";
      parameter Medium.SpecificEnthalpy h_start 
        "Start value of specific enthalpy (used as enthalpy for topPorts if back flow)";
      parameter Medium.MassFraction X_start[Medium.nX] 
        "Start value of mass fractions m_i/m";
      parameter Modelica.SIunits.Height level_start "Start value of tank level";
      
      parameter Modelica.SIunits.AbsolutePressure p_ambient 
        "Tank surface pressure"                                   annotation(Dialog);
      input Modelica.SIunits.Height level "Actual tank level" 
                                                annotation(Dialog);
      input Medium.SpecificEnthalpy h 
        "Actual specific enthalpy of fluid in tank"               annotation(Dialog);
      input Medium.Density d "Actual specific density of fluid in tank" 
                                                        annotation(Dialog);
      input Medium.MassFraction Xi[Medium.nXi] 
        "Actual mass fractions of fluid in tank"                  annotation(Dialog);
      parameter Real k_small(min=0) = 1e-5 
        "Small regularization range if tank level is below bottom_height or side_height; k_small = 0 gives ideal switch"
                annotation(Evaluate=true);
      parameter Real s_start = 0;
      
      output Medium.EnthalpyFlowRate H_flow 
        "= port.H_flow (used to transform vector of connectors in vector of Real numbers)";
      output Medium.MassFlowRate m_flow 
        "= port.m_flow (used to transform vector of connectors in vector of Real numbers)";
      output Medium.MassFlowRate mXi_flow[Medium.nXi] 
        "= port.mXi_flow (used to transform vector of connectors in vector of Real numbers)";
      
    annotation (Documentation(info="<html>
<p>
This component contains the equations that attach the pipe
to the tank. The main reason to introduce this component is
that Dymola has currently limitations for connector arrays
when the dimension is zero. Without this utility component
it would not be possible to set, e.g., n_topPorts to zero.
</p>
 
 
</html>"), Icon(Rectangle(extent=[-100,0; 100,-100], style(
            color=3,
            rgbcolor={0,0,255},
            fillColor=7,
            rgbfillColor={255,255,255})), Text(
          extent=[-122,48; 132,6],
          style(
            color=3,
            rgbcolor={0,0,255},
            fillColor=7,
            rgbfillColor={255,255,255},
            fillPattern=1),
          string="%name")));
    protected 
    outer Modelica_Fluid.Ambient ambient "Global default options";
    parameter Modelica.SIunits.Area pipeArea=Modelica.Constants.pi*(
          pipeDiameter/2)^2;
    parameter Medium.MassFlowRate m_flow_nominal = 1 
        "Nominal mass flow rate used for scaling (has only an effect if k_small > 0)";
    parameter Medium.AbsolutePressure p_nominal = p_ambient;
    Modelica.SIunits.Length aboveLevel=level - pipeHeight;
    Boolean m_flow_out(start=true,fixed=true) "true= massflow out of tank";
    Real s(start=s_start) 
        "path parameter of parameterized curve description (either m_flow/m_flow_nominal or (port.p-p_ambient)/p_ambient)";
  equation 
    H_flow = port.H_flow;
    m_flow = port.m_flow;
    mXi_flow = port.mXi_flow;
      
    if onlyInFlow then
       m_flow_out = false "Dummy value in this if branch";
       port.p = p_ambient;
       /* flow should never out of the port. However, if this occurs in a 
        small time interval (e.g. during initialization), the start values of
        h and X are provided, since otherwise there is a singular
        system 
     */
       port.H_flow = semiLinear(port.m_flow, port.h, h_start);
       port.mXi_flow = semiLinear(port.m_flow, port.Xi, X_start[1:Medium.nXi]);
       assert(port.m_flow > -1e-6, "Mass flows out of tank via topPort. This indicates a wrong model");
       s = 0;
    else
       port.H_flow = semiLinear(port.m_flow, port.h, h);
       port.mXi_flow = semiLinear(port.m_flow, port.Xi, Xi);
        
  /* Original equations from Poschlad/Remelhe:
*/
       s = 0;
       m_flow_out = (pre(m_flow_out) and not port.p>p_ambient) or (port.m_flow < -1e-6);
        
       if (aboveLevel > 0) then
         port.p = aboveLevel*ambient.g*d + p_ambient - smooth(2,noEvent(if m_flow < 0 then m_flow^2/(2*d*pipeArea^2) else 0));
       else
         if pre(m_flow_out) then
            m_flow = 0;
         else
            port.p = p_ambient;
         end if;
       end if;
        
  /* Martin Otter: The following equations are a declarative form 
   (parameterized curve description) of the above equations and
   should theoretically work better. However, some examples with
   IF97 water fail, whereas the above works. Therefore, not used. 
       Add the following text to OpenTank, once the initialization
   with this solution works for Modelica_Fluid.Examples.Tanks.ThreeOpenTanks:
 
OpenTank:
<p>
The situation when the tank level is below bottom_heights[i] or side_heights[i]
is handeled properly. Details are described
<a href="Modelica:Modelica_Fluid.Utilities.TankAttachment">here</a>
</p> 
 
TankAttachment:
<p>
If a bottom or side connector is above the actual tank level, the
following characteristic is used to compute the mass flow rate port.m_flow
from the connector to the tank and the absolute pressure port.p
in the port:
</p>
 
<img src="../Images/Components/Tank_PipeAboveTankLevel.png">   
 
 
     m_flow_out = s <= 0;
     if aboveLevel >= 0 then
        m_flow = m_flow_nominal*s "equation to compute s, which is a dummy in this branch";
        port.p - p_ambient = aboveLevel*ambient.g*d  -
                             smooth(2,if m_flow_out then s*abs(s)*m_flow_nominal^2/(2*d*pipeArea^2) else k_small*m_flow_nominal*s);
     else
        m_flow = (if m_flow_out then k_small*p_nominal else m_flow_nominal)*s;
        port.p - p_ambient = (if m_flow_out then p_nominal else k_small*m_flow_nominal)*s;
     end if;
*/
        
    end if;
      
    /*
  More precise equations (introduce them later; need to transform
  them from energy balance form 1 to form 2):
 
  Momentum balance:
  (integrated momentum equation for frictionless fluid with density that is
   independent of the level, i.e., the unsteady Bernoulli equation for incompressible fluid)
  v_level = der(level);
  v = -port.m_flow/(rho*A_outlet);
  level*der(v_level) + (v^2 - v_level^2)/2 - g*level + (p - p_ambient)/rho = 0; or
  rho*level*der(v_level) + rho*(v^2 - v_level^2)/2 - rho*g*level + (p - p_ambient) = 0;
 
  Energy balance:
  Potential energy: E_pot = integ(dm*g*s)
                          = g*integ(rho*A*s*ds)
                          = g*rho*A*z^2/2
  Kinetic energy  : E_kin = integ(dm*v^2/2)
                          = integ(rho*A*v^2/2*ds)
                          = rho*A*v^2/2*integ(ds)
                          = rho*A*v^2/2*z
                          = M*v^2/2
  E = U + M*g*z/2 + M*v_level^2/2
  der(E) = port.H_flow + port.m_flow*v^2/2 - p_ambient*area*der(level)
*/
      
  end TankAttachment;
    
   record PressureLossFactors 
      "Data structure defining constant loss factors zeta for dp = zeta*rho*v*|v|/2 and functions providing the data for some loss types" 
      extends Modelica.Icons.Record;
      
    parameter Modelica.SIunits.Diameter D_a "Diameter at port_a";
    parameter Modelica.SIunits.Diameter D_b=D_a "Diameter at port_b";
    parameter Real zeta1 "Loss factor for flow port_a -> port_b";
    parameter Real zeta2=zeta1 "Loss factor for flow port_b -> port_a";
    parameter Modelica.SIunits.ReynoldsNumber Re_turbulent=1000 
        "Loss factors suited for |Re| >= Re_turbulent";
    parameter Modelica.SIunits.Diameter D_Re=min(D_a, D_b) 
        "Diameter used to compute Re";
    parameter Boolean zeta1_at_a = true 
        "dp = zeta1*(if zeta1_at_a then d_a*v_a^2/2 else d_b*v_b^2/2)";
    parameter Boolean zeta2_at_a = false 
        "dp = -zeta2*(if zeta2_at_a then d_a*v_a^2/2 else d_b*v_b^2/2)";
    parameter Boolean zetaLaminarKnown = false 
        "= true, if zeta = c0/Re in laminar region";
    parameter Real c0 = 1 
        "zeta = c0/Re; dp = zeta*d_Re*v_Re^2/2, Re=v_Re*D_Re*d_Re/eta_Re)"                   annotation(Dialog(enable=zetaLaminarKnown));
      
    annotation (preferedView="info", Documentation(info="<html>
<p>
This record defines the pressure loss factors of a pipe
segment (orifice, bending etc.) with a minimum amount of data.
If available, data should be provided for <b>both flow directions</b>,
i.e., flow from port_a to port_b and from port_b to port_a, 
as well as for the <b>laminar</b> and the <b>turbulent</b> region.
It is also an option to provide the loss factor <b>only</b> for the
<b>turbulent</b> region for a flow from port_a to port_b.
</p>
<p>
The following equations are used:
</p>
<pre>   &Delta;p = 0.5*&zeta;*&rho;*v*|v|
      = 0.5*&zeta;/A^2 * (1/&rho;) * m_flow*|m_flow|
        Re = |v|*D*&rho;/&eta;
</pre>
<table border=1 cellspacing=0 cellpadding=2>
<tr><td><b>flow type</b></td>
    <td><b>&zeta;</b> = </td>
    <td><b>flow region</b></td></tr>
<tr><td>turbulent</td>
    <td><b>zeta1</b> = const.</td>
    <td>Re &ge;  Re_turbulent, v &ge; 0</td></tr>
<tr><td></td>
    <td><b>zeta2</b> = const.</td>
    <td>Re &ge; Re_turbulent, v &lt; 0</td></tr>
<tr><td>laminar</td>
    <td><b>c0</b>/Re</td>
    <td>both flow directions, Re small; c0 = const.</td></tr>
</table>
<p>
where
</p>
<ul>
<li> &Delta;p is the pressure drop: &Delta;p = port_a.p - port_b.p</li>
<li> v is the mean velocity.</li>
<li> &rho; is the density.</li>
<li> &zeta; is the loss factor that depends on the geometry of
     the pipe. In the turbulent flow regime, it is assumed that
     &zeta; is constant and is given by \"zeta1\" and
     \"zeta2\" depending on the flow direction.<br>
     When the Reynolds number Re is below \"Re_turbulent\", the
     flow is laminar for small flow velocities. For higher 
     velocities there is a transition region from 
     laminar to turbulent flow. The loss factor for
     laminar flow at small velocities is defined by the often occuring
     approximation c0/Re. If c0 is different for the two
     flow directions, the mean value has to be used 
     (c0 = (c0_ab + c0_ba)/2).<li>
<li> The equation \"&Delta;p = 0.5*&zeta;*&rho;*v*|v|\" is either with
     respect to port_a or to port_b, depending on the definition
     of the particular loss factor &zeta; (in some references loss
     factors are defined with respect to port_a, in other references
     with respect to port_b).</li>
 
<li> Re = |v|*D_Re*&rho;/&eta; = |m_flow|*D_Re/(A_Re*&eta;) 
     is the Reynolds number at the smallest cross
     section area. This is often at port_a or at port_b, but can
     also be between the two ports. In the record, the diameter
     D_Re of this smallest cross section area has to be provided, as
     well, as Re_turbulent, the absolute value of the 
     Reynolds number at which
     the turbulent flow starts. If Re_turbulent is different for
     the two flow directions, use the smaller value as Re_turbulent.</li>
<li> D is the diameter of the pipe. If the pipe has not a 
     circular cross section, D = 4*A/P, where A is the cross section
     area and P is the wetted perimeter.</li>
<li> A is the cross section area with A = &pi;(D/2)^2.
<li> &eta; is the dynamic viscosity.</li>
</ul>
<p>
The laminar and the transition region is usually of
not much technical interest because the operating point is
mostly in the turbulent regime. For simplification and for
numercial reasons, this whole region is described by two
polynomials of third order, one polynomial for m_flow &ge; 0 
and one for m_flow &lt; 0. The polynomials start at 
Re = |m_flow|*4/(&pi;*D_Re*&eta;), where D_Re is the
smallest diameter between port_a and port_b.
The common derivative
of the two polynomials at Re = 0 is
computed from the equation \"c0/Re\". Note, the pressure drop
equation above in the laminar region is always defined
with respect to the smallest diameter D_Re.
</p>
<p>
If no data for c0 is available, the derivative at Re = 0 is computed in such
a way, that the second derivatives of the two polynomials
are identical at Re = 0. The polynomials are constructed, such that
they smoothly touch the characteristic curves in the turbulent
regions. The whole characteristic is therefore <b>continuous</b>
and has a <b>finite</b>, <b>continuous first derivative everywhere</b>.
In some cases, the constructed polynomials would \"vibrate\". This is 
avoided by reducing the derivative at Re=0 in such a way that
the polynomials are guaranteed to be monotonically increasing.
The used sufficient criteria for monotonicity follows from:
</p>
 
<dl>
<dt> Fritsch F.N. and Carlson R.E. (1980):</dt>
<dd> <b>Monotone piecewise cubic interpolation</b>.
     SIAM J. Numerc. Anal., Vol. 17, No. 2, April 1980, pp. 238-246</dd>
</dl>
</html>"));
      
     encapsulated function wallFriction 
        "Compute pressure loss factors due to friction in a straight pipe with walls of nonuniform roughness (commercial pipes)" 
        import Modelica_Fluid_WorkInProgress.Utilities.PressureLossFactors;
        import lg = Modelica.Math.log10;
        import SI = Modelica.SIunits;
        
       input SI.Length length "Length of pipe";
       input SI.Diameter diameter "Inner diameter of pipe";
       input SI.Length roughness(min=1e-10) 
          "Absolute roughness of pipe (> 0 required, details see info layer)";
       output PressureLossFactors data 
          "Pressure loss factors for both flow directions";
       annotation (Icon(Rectangle(extent=[-100,48; 100,-50], style(
               color=0,
               rgbcolor={0,0,0},
               fillColor=7,
               rgbfillColor={255,255,255}))),
                                 Diagram(
           Rectangle(extent=[-100,48; 100,-50], style(
                 color=0,
                 rgbcolor={0,0,0},
                 fillColor=7,
                 rgbfillColor={255,255,255},
                 fillPattern=1)),
           Line(points=[-60,-50; -60,48], style(
               color=3,
               rgbcolor={0,0,255},
               arrow=3,
               fillColor=3,
               rgbfillColor={0,0,255},
               fillPattern=1)),
           Text(
             extent=[-50,16; 6,-10],
             style(
               color=3,
               rgbcolor={0,0,255},
               arrow=3,
               fillColor=3,
               rgbfillColor={0,0,255},
               fillPattern=1),
             string="diameter"),
           Line(points=[-100,60; 100,60], style(
               color=3,
               rgbcolor={0,0,255},
               arrow=3,
               fillColor=3,
               rgbfillColor={0,0,255},
               fillPattern=1)),
           Text(
             extent=[-34,80; 34,62],
             style(
               color=3,
               rgbcolor={0,0,255},
               arrow=3,
               fillColor=3,
               rgbfillColor={0,0,255},
               fillPattern=1),
             string="length")),
         Documentation(info="<html>
<p>
Friction in straight pipe with walls of nonuniform roughness 
(commercial pipes)
</p>
<p>
The loss factors are given for mass flow rates from 
port_a to port_b as:
</p>
<pre>
  turbulent flow (Idelchik 1994, diagram 2-5, p. 117)
     zeta = (L/D)/(2*lg(3.7 / &Delta;))^2, for Re >= 560/&Delta;
     for Re >= 560/&Delta; the loss factor does not depend on the
     Reynolds number. For Re >= 4000, the flow is turbulent,
     but depends both on &Delta; and slightly on Re.
&nbsp;
  laminar flow (Idelchick 1994, diagram 2-1, p. 110):
     zeta = 64*(L/D)/Re
</pre>
<p>
where
</p>
<ul>
<li> D is the inner pipe diameter</li>
<li> L is the lenght of the pipe</li>
<li> &Delta; = &delta;/D is the relative roughness where &delta; is
     the absolute \"roughness\", i.e., the averaged height of asperities in the pipe.
     (&delta; may change over time due to growth of surface asperities during
      service, see [Idelchick 1994, p. 85, Tables 2-1, 2-2]).</li>
</ul>
<p>
The absolute roughness <font face=\"Symbol\">d</font> has usually to
be estimated. In <i>[Idelchik 1994, pp. 105-109,
Table 2-5; Miller 1990, p. 190, Table 8-1]</i> many examples are given.
As a short summary:
</p>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><td><b>Smooth pipes</b></td>
      <td>Drawn brass, coper, aluminium, glass, etc.</td>
      <td><font face=\"Symbol\">d</font> = 0.0025 mm</td>
  </tr>
  <tr><td rowspan=\"3\"><b>Steel pipes</b></td>
      <td>New smooth pipes</td>
      <td><font face=\"Symbol\">d</font> = 0.025 mm</td>
  </tr>
  <tr><td>Mortar lined, average finish</td>
      <td><font face=\"Symbol\">d</font> = 0.1 mm</td>
  </tr>
  <tr><td>Heavy rust</td>
      <td><font face=\"Symbol\">d</font> = 1 mm</td>
  </tr>
  <tr><td rowspan=\"3\"><b>Concrete pipes</b></td>
      <td>Steel forms, first class workmanship</td>
      <td><font face=\"Symbol\">d</font> = 0.025 mm</td>
  </tr>
  <tr><td>Steel forms, average workmanship</td>
      <td><font face=\"Symbol\">d</font> = 0.1 mm</td>
  </tr>
  <tr><td>Block linings</td>
      <td><font face=\"Symbol\">d</font> = 1 mm</td>
  </tr>
</table>
</html>"));
      protected 
       Real Delta = roughness/diameter "relative roughness";
     algorithm 
       data.D_a          := diameter;
       data.D_b          := diameter;
       data.zeta1        := (length/diameter)/(2*lg(3.7 /Delta))^2;
       data.zeta2        := data.zeta1;
       data.Re_turbulent := 4000 
          ">= 560/Delta flow does not depend on Re, but interpolation is bad";
       data.D_Re         := diameter;
       data.zeta1_at_a   := true;
       data.zeta2_at_a   := false;
       data.zetaLaminarKnown := true;
       data.c0               := 64*(length/diameter);
     end wallFriction;
      
     encapsulated function suddenExpansion 
        "Compute pressure loss factors for sudden expansion or contraction in a pipe (for both flow directions)" 
        import Modelica_Fluid_WorkInProgress.Utilities.PressureLossFactors;
        import SI = Modelica.SIunits;
        
       input SI.Diameter D_a "Inner diameter of pipe at port_a";
       input SI.Diameter D_b "Inner diameter of pipe at port_b";
       output PressureLossFactors data 
          "Pressure loss factors for both flow directions";
       annotation (Icon(
           Rectangle(extent=[-100,40; 0,-40], style(
               color=7,
               rgbcolor={255,255,255},
               fillColor=7,
               rgbfillColor={255,255,255})),
           Rectangle(extent=[0,100; 100,-100], style(
               color=7,
               rgbcolor={255,255,255},
               fillColor=7,
               rgbfillColor={255,255,255})),
           Line(points=[0,40; -100,40; -100,-40; 0,-40; 0,-100; 100,-100; 100,100; 0,
                 100; 0,40], style(
               color=0,
               rgbcolor={0,0,0},
               fillColor=7,
               rgbfillColor={255,255,255},
               fillPattern=1))), Diagram(
           Line(points=[0,40; -100,40; -100,-40; 0,-40; 0,-100; 100,-100; 100,100; 0,
                 100; 0,40], style(
               color=0,
               rgbcolor={0,0,0},
               fillColor=7,
               rgbfillColor={255,255,255},
               fillPattern=1)),
           Rectangle(extent=[-100,40; 0,-40], style(
               color=7,
               rgbcolor={255,255,255},
               fillColor=7,
               rgbfillColor={255,255,255})),
           Rectangle(extent=[0,100; 100,-100], style(
               color=7,
               rgbcolor={255,255,255},
               fillColor=7,
               rgbfillColor={255,255,255})),
           Line(points=[0,40; -100,40; -100,-40; 0,-40; 0,-100; 100,-100; 100,100; 0,
                 100; 0,40], style(
               color=0,
               rgbcolor={0,0,0},
               fillColor=7,
               rgbfillColor={255,255,255},
               fillPattern=1)),
           Line(points=[-60,-40; -60,40], style(
               color=3,
               rgbcolor={0,0,255},
               arrow=3,
               fillColor=3,
               rgbfillColor={0,0,255},
               fillPattern=1)),
           Text(
             extent=[-50,16; -26,-10],
             style(
               color=3,
               rgbcolor={0,0,255},
               arrow=3,
               fillColor=3,
               rgbfillColor={0,0,255},
               fillPattern=1),
             string="D_a"),
           Line(points=[34,-100; 34,100], style(
               color=3,
               rgbcolor={0,0,255},
               arrow=3,
               fillColor=3,
               rgbfillColor={0,0,255},
               fillPattern=1)),
           Text(
             extent=[54,16; 78,-10],
             style(
               color=3,
               rgbcolor={0,0,255},
               arrow=3,
               fillColor=3,
               rgbfillColor={0,0,255},
               fillPattern=1),
             string="D_b")),
         Documentation(info="<html>
<p>
The loss factors are given for mass flow rates from 
port_a to port_b as:
</p>
<pre>
   A_a &lt; A_b (Idelchik 1994, diagram 4-1, p. 208):
      zeta = dp/(d_a*v_a^2/2)
           = (1 - A_a/A_b)^2 for Re_a &ge; 3.3e3 (turbulent flow)
      zeta = 30/Re           for Re_a &lt; 10    (laminar flow)
&nbsp;
   A_a &gt; A_b (Idelchik 1994, diagram 4-9, p. 216 and diagram 4-10, p. 217)
      zeta = dp/(d_b*v_b^2/2)
           = 0.5*(1 - A_b/A_a)^0.75 for Re_b &ge; 1e4 (turbulent flow)
      zeta = 30/Re                  for Re_a &lt; 10  (laminar flow)
</pre>
</html>"));
      protected 
       Real A_rel;
     algorithm 
       data.D_a          := D_a;
       data.D_b          := D_b;
       data.Re_turbulent := 100;
       data.zetaLaminarKnown := true;
       data.c0 := 30;
        
       if D_a <= D_b then
          A_rel :=(D_a/D_b)^2;
          data.zeta1 :=(1 - A_rel)^2;
          data.zeta2 :=0.5*(1 - A_rel)^0.75;
          data.zeta1_at_a :=true;
          data.zeta2_at_a :=true;
          data.D_Re := D_a;
       else
          A_rel :=(D_b/D_a)^2;
          data.zeta1 :=0.5*(1 - A_rel)^0.75;
          data.zeta2 :=(1 - A_rel)^2;
          data.zeta1_at_a :=false;
          data.zeta2_at_a :=false;
          data.D_Re := D_b;
       end if;
     end suddenExpansion;
      
     encapsulated function sharpEdgedOrifice 
        "Compute pressure loss factors for sharp edged orifice (for both flow directions)" 
        import NonSI = Modelica.SIunits.Conversions.NonSIunits;
        import Modelica_Fluid_WorkInProgress.Utilities.PressureLossFactors;
        import SI = Modelica.SIunits;
        
       input SI.Diameter D_pipe 
          "Inner diameter of pipe (= same at port_a and port_b)";
       input SI.Diameter D_min "Smallest diameter of orifice";
       input SI.Diameter L "Length of orifice";
       input NonSI.Angle_deg alpha "Angle of orifice";
       output PressureLossFactors data 
          "Pressure loss factors for both flow directions";
       annotation (Icon(Rectangle(extent=[-100,60; 100,-60], style(
                 color=0,
                 rgbcolor={0,0,0},
                 fillColor=7,
                 rgbfillColor={255,255,255})),
             Polygon(points=[-30,60; -30,12; 30,50; 30,60; -30,60], style(
                 color=0,
                 rgbcolor={0,0,0},
                 fillColor=7,
                 rgbfillColor={255,255,255},
                 fillPattern=8)),
             Polygon(points=[-30,-10; -30,-60; 30,-60; 30,-50; -30,-10], style(
                 color=0,
                 rgbcolor={0,0,0},
                 fillColor=7,
                 rgbfillColor={255,255,255},
                 fillPattern=8))),
                                 Diagram(
                        Rectangle(extent=[-100,60; 100,-60], style(
                 color=0,
                 rgbcolor={0,0,0},
                 fillColor=7,
                 rgbfillColor={255,255,255})),
             Polygon(points=[-30,60; -30,12; 30,50; 30,60; -30,60], style(
                 color=0,
                 rgbcolor={0,0,0},
                 fillColor=7,
                 rgbfillColor={255,255,255},
                 fillPattern=8)),
             Polygon(points=[-30,-10; -30,-60; 30,-60; 30,-50; -30,-10], style(
                 color=0,
                 rgbcolor={0,0,0},
                 fillColor=7,
                 rgbfillColor={255,255,255},
                 fillPattern=8)),
           Line(points=[-82,-60; -82,60], style(
               color=3,
               rgbcolor={0,0,255},
               arrow=3,
               fillColor=3,
               rgbfillColor={0,0,255},
               fillPattern=1)),
           Text(
             extent=[-78,16; -44,-8],
             style(
               color=3,
               rgbcolor={0,0,255},
               arrow=3,
               fillColor=3,
               rgbfillColor={0,0,255},
               fillPattern=1),
               string="D_pipe"),
           Line(points=[-30,-10; -30,12], style(
               color=3,
               rgbcolor={0,0,255},
               arrow=3,
               fillColor=3,
               rgbfillColor={0,0,255},
               fillPattern=1)),
           Text(
             extent=[-24,14; 8,-10],
             style(
               color=3,
               rgbcolor={0,0,255},
               arrow=3,
               fillColor=3,
               rgbfillColor={0,0,255},
               fillPattern=1),
               string="D_min"),
           Text(
             extent=[-20,84; 18,70],
             style(
               color=3,
               rgbcolor={0,0,255},
               arrow=3,
               fillColor=3,
               rgbfillColor={0,0,255},
               fillPattern=1),
               string="L"),
           Line(points=[30,68; -30,68],   style(
               color=3,
               rgbcolor={0,0,255},
               arrow=3,
               fillColor=3,
               rgbfillColor={0,0,255},
               fillPattern=1)),
             Line(points=[16,40; 32,18; 36,-2; 34,-20; 20,-42], style(
                 color=3,
                 rgbcolor={0,0,255},
                 arrow=3,
                 fillColor=3,
                 rgbfillColor={0,0,255},
                 fillPattern=8)),
             Text(
               extent=[38,8; 92,-6],
               style(
                 color=3,
                 rgbcolor={0,0,255},
                 arrow=3,
                 fillColor=3,
                 rgbfillColor={0,0,255},
                 fillPattern=8),
               string="alpha")),
         Documentation(info="<html>
<p>
Loss factor for mass flow rate from port_a to port_b
(Idelchik 1994, diagram 4-14, p. 221):
</p>
<pre>
   zeta = [(1-A0/A1) + 0.707*(1-A0/A1)^0.375]^2*(A1/A0)^2 
          for Re(A0) >= 1e5,  independent of alpha
</pre>
<p>
Loss factor for mass flow rate from port_b to port_a
(Idelchik 1994, diagram 4-13, p. 220, with A2=A1):
</p>
<pre>
   zeta = k*(1 - A0/A1)^0.75 + (1 - A0/A1)^2 + 2*sqrt(k*(1-A0/A1)^0.375) + (1- A0/A1)
          k  = 0.13 + 0.34*10^(-(3.4*LD+88.4*LD^2.3)) 
               (there is a typing error in the formula in diagram 4-13, the above
                equation corresponds to table (a) in diagram 4-12)
          LD = L/D0
          for Re(A0) >= 1e4, 40 deg &le; alpha &le; 60 deg 
                             for other values of alpha, k is given as table
                             in diagram 3-7 (this is not yet included in the function)
</pre
</html>"));
      protected 
       Real D_rel = D_min/D_pipe;
       Real LD = L/D_min;
       Real k = 0.13 + 0.34*10^(-(3.4*LD+88.4*LD^2.3));
     algorithm 
       data.D_a          := D_pipe;
       data.D_b          := D_pipe;
       data.zeta1        := ((1-D_rel) + 0.707*(1-D_rel)^0.375)^2*(1/D_rel)^2;
       data.zeta2        := k*(1 - D_rel)^0.75 + (1 - D_rel)^2 +
                            2*sqrt(k*(1-D_rel)^0.375) + (1- D_rel);
       data.Re_turbulent := 1e4;
       data.D_Re         := D_min;
       data.zeta1_at_a   := true;
       data.zeta2_at_a   := false;
       data.zetaLaminarKnown := false;
       data.c0               := 0;
     end sharpEdgedOrifice;
   end PressureLossFactors;
    
    function diameter_of_squarePipe 
      "Determine hydraulic diameter of pipe with square cross sectional area" 
        extends Modelica.Icons.Function;
      input Modelica.SIunits.Length width "Inner width of pipe";
      input Modelica.SIunits.Length height " Inner height of pipe";
      output Modelica.SIunits.Diameter D "Inner (hydraulic) diameter of pipe";
    algorithm 
      D :=4*width*height/(2*(width+height));
    end diameter_of_squarePipe;
    
      model FiniteVolume 
      "One dimensional volume according to the finite volume method with 1 mass, 1 energy and 2 momentum balances" 
      
      import Modelica.Math;
      
        replaceable package Medium = 
          Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
            annotation (choicesAllMatching = true);
      
        Modelica_Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium 
          =                                                                     Medium) 
          annotation(extent=[-120, -10; -100, 10]);
        Modelica_Fluid.Interfaces.FluidPort_b port_b(redeclare package Medium 
          =                                                                      Medium) 
          annotation(extent=[120, -10; 100, 10]);
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort 
          annotation(extent=[-10, 60; 10, 80], rotation=-90);
        Medium.BaseProperties medium(preferredMediumStates=true) 
        "Medium properties in the middle of the finite volume";
        Modelica.SIunits.Mass M "Total mass in volume";
        Modelica.SIunits.Mass[Medium.nXi] MXi "Independent component masses";
        Modelica.SIunits.Energy U "Inner energy";
        parameter Modelica.SIunits.Length L "Length of volume";
      
        parameter Modelica.SIunits.Area A_a;
        parameter Modelica.SIunits.Area A_b=A_a;
      
        parameter Modelica.SIunits.Length Z_a=0;
        parameter Modelica.SIunits.Length Z_b=Z_a;
      
        parameter Boolean dynamicMomentumBalance=false 
        "If false, der(m_flow) is neglected in momentum balance" 
                                                       annotation(Evaluate=true,
            Dialog(tab="Level of Detail"));
        parameter Boolean includeKineticTerm=false 
        "If false, d*v^2 is neglected in momentum balance" 
                                                   annotation(Evaluate=true,
            Dialog(tab="Level of Detail"));
        parameter Boolean includeViscosity=false 
        "If false, artifical viscosity is neglected" 
                                                annotation(Evaluate=true, Dialog(tab=
                "Level of Detail"));
        parameter Real viscosityFactor1=0 annotation(Dialog(enable=includeViscosity,tab="Level of Detail"));
        parameter Real viscosityFactor2=1 annotation(Dialog(enable=includeViscosity,tab="Level of Detail"));
      
        input Modelica.SIunits.Pressure dp_a 
        "Pressure loss due to pipe friction between port_a and middle of pipe";
        input Modelica.SIunits.Pressure dp_b 
        "Pressure loss due to pipe friction between middle of pipe and port_b";
      
        Real my "Artifical viscosity";
      
        annotation (
          Documentation(info="<html>
<p>
Model <b>FiniteVolume</b> is a generic finite volume for 1-dim. thermo-fluid flow
in piping networks which has the following properties:
</p>
<ul>
<li> A FiniteVolume model is <b>independent</b> of the <b>medium</b> model. 
     The only requirement is that the medium model has to be
     a subclass of Modelica.Media.Interfaces.<b>PartialMedium</b>.
     As a consequence, the FiniteVolume model can be used
     for incompressible or compressible media, fluids with one 
     and multiple substances as well as for one and multiple phases. 
     The FiniteVolume model depends also not on the independent
     variables of the medium model.</li>
<li> A FiniteVolume model contains a <b>staggered grid</b> consisting
     of one volume from port_a to port_b for which the mass and
     energy balances are formulated and from two volumes from port_a 
     to the middle of the FiniteVolume and from the middle to
     port_b for which momentum
     balances are provided. When a FiniteVolume model is connected
     to another FiniteVolume, the adjacent momentum balance volumes
     are automatically merged together to form one volume. Thus, 
     a connected network of FiniteVolumes consists of a staggered
     grid for the balance equations.</li>
<li> For the intensive properties, such as density, temperature,
     an <b>upwind scheme</b> is used, i.e., depending on the direction
     of the mass flow rate, the upwind value is used in the equations.
     Also zero mass flow rate is handeled appropriately.</li>
<li> In order that the FiniteVolume model can be utilized, equations
     for the pipe friction have to be added via the input variables
     dp_a (pressure loss from port_a to the middle of the FiniteVolume)
     and dp_b (pressure loss from port_b to the middle of the FiniteVolume).</li>
<li> A FiniteVolume component contains <b>one medium</b> model in the 
     <b>middle</b> of the FiniteVolume.</li>
</ul>
</html>",       revisions="<html>
<ul>
<li><i>Aug. 30, 2004</i>
    by Hilding Elmqvist, Dynasim:<br>
    Further improvements + artifical viscosity introduced to remove
    unphysical oscillations in shock waves.</li>
<li><i>May 28, 2004</i>
    by Hilding Elmqvist, Dynasim:<br>
    Implemented.</li>
</ul>
</html>"),Diagram,
          Icon(Rectangle(extent=[-100, -60; 100, 60], style(
                color=3,
                rgbcolor={0,0,255},
                gradient=2,
                fillColor=76,
                rgbfillColor={170,170,255})),
            Ellipse(extent=[-16,16; 14,-12],    style(fillColor=0)),
            Rectangle(extent=[-90,46; 92,-46], style(color=0, rgbcolor={0,0,0})),
            Rectangle(extent=[-80,36; -6,-36], style(color=0, rgbcolor={0,0,0})),
            Rectangle(extent=[8,36; 82,-36], style(color=0, rgbcolor={0,0,0}))));
    protected 
        Modelica.SIunits.MassFlowRate m_flow_a;
        Modelica.SIunits.MassFlowRate m_flow_b;
        Modelica.SIunits.MassFlowRate m_flow_middle;
        constant Real pi=Modelica.Constants.pi;
        constant Real g=Modelica.Constants.g_n;
        parameter Modelica.SIunits.Area A_m=(A_a + A_b)/2;
        parameter Modelica.SIunits.Length dx=L;
      equation 
        //Extensive properties
          M=medium.d*A_m*dx;
          MXi=M*medium.Xi;
          U=M*medium.u;
      
        // Mass balance over the interval a to b
        //der(medium.d)*A_m*dx = port_a.m_flow + port_b.m_flow;
        der(M)=port_a.m_flow + port_b.m_flow;
      
        // Substance mass balances over the interval a to b
        // der(medium.d*medium.X)*A_m*dx = port_a.mXi_flow + port_b.mXi_flow;
        //(der(medium.d)*medium.X + medium.d*der(medium.X))*A_m*dx = port_a.mXi_flow + port_b.mXi_flow;
        der(MXi)= port_a.mXi_flow + port_b.mXi_flow;
      
        // Energy balance over the interval a to b
        // der(medium.d*medium.u)*A_m*dx = port_a.H_flow + port_b.H_flow + m_flow_middle/
        //   medium.d*(port_b.p - port_a.p) + heatPort.Q_flow;
        //(der(medium.d)*medium.u + medium.d*der(medium.u))*A_m*dx = port_a.H_flow + port_b.H_flow - m_flow_middle/
        //  medium.d*(port_a.p - port_b.p - dp_a - dp_b) + heatPort.Q_flow;
        der(U)= port_a.H_flow + port_b.H_flow - m_flow_middle/medium.d*(port_a.p - port_b.p - dp_a - dp_b) + heatPort.Q_flow;
      
        m_flow_middle = (port_a.m_flow - port_b.m_flow)/2 
        "since assumed same density in entire interval a to b";
      
        // Momentum balance over interval a to dx/2
        (if dynamicMomentumBalance then der(m_flow_a)*dx/2 else 0) =
          A_m*(port_a.p - medium.p - dp_a) +
          (if includeKineticTerm then 
            - m_flow_middle^2/(A_m*medium.d) else 0)
          - A_m*medium.d/2*g*(Z_b - Z_a) +
          (if includeViscosity then my*((-m_flow_b-m_flow_a)/dx - 0) else 0);
          /* Removed: port_a.m_flow^2/(A_a*medium_a.d) */
      
        // Momentum balance over interval dx/2 to b
        (if dynamicMomentumBalance then -der(m_flow_b)*dx/2 else 0) =
          A_m*(medium.p - port_b.p - dp_b) +
          (if includeKineticTerm then 
            m_flow_middle^2/(A_m*medium.d) else 0)
          - A_m*medium.d/2*g*(Z_b - Z_a) +
          (if includeViscosity then my*(0 - (-m_flow_b-m_flow_a)/dx) else 0);
          /* Removed: - port_b.m_flow^2/(A_b*medium_b.d) */
      
         if includeViscosity then
           my = viscosityFactor1 + viscosityFactor2*dx*(if m_flow_middle*(-m_flow_b-m_flow_a) < 0 then 
              abs(-m_flow_b-m_flow_a)/(A_m*medium.d) else 0);
         else
           my = 0;
         end if;
      
        // Coupling to environment  
        m_flow_a = port_a.m_flow 
        "Due to problem with non-aliasing and semiLinear";
        m_flow_b = port_b.m_flow;
      
        // Upwind scheme (use properties from upwind port and handle zero flow)  
        port_a.H_flow = semiLinear(port_a.m_flow, port_a.h, medium.h);
        port_b.H_flow = semiLinear(port_b.m_flow, port_b.h, medium.h);
        port_a.mXi_flow = semiLinear(port_a.m_flow, port_a.Xi, medium.Xi);
        port_b.mXi_flow = semiLinear(port_b.m_flow, port_b.Xi, medium.Xi);
      
        // Heat port has the medium temperature
        heatPort.T = medium.T;
      
      end FiniteVolume;
    
  model PipeSegment 
      "One segment of a pipe with 1 mass, 1 energy, 2 momementum balances and pipe friction" 
      
      import Modelica_Fluid.Types.Init.*;
    extends FiniteVolume(medium(
               p(start=p_start),
               T(start=T_start), h(start=h_start), Xi(start=X_start[1:Medium.nXi])));
    extends 
        Modelica_Fluid_WorkInProgress.Interfaces.PartialInitializationParameters;
      
    parameter Boolean linearPressureDrop=true;
    parameter Modelica.SIunits.AbsolutePressure dp_nominal(min=1.e-10)=1 
        "Nominal pressure drop";
    parameter Modelica.SIunits.MassFlowRate m_flow_nominal=1E-3 
        "Nominal mass flow rate at nominal pressure drop";
      
    annotation (Documentation(info="<html>
<p>
Model <b>PipeSegment</b> describes one segment of a pipe.
It consists of the following parts:
</p>
<ul>
<li> One <a href=\"Modelica:Modelica_Fluid.Utilities.FiniteVolume\">FiniteVolume</a>
     model described by 1 mass, 1 energy, and 2 momemtum balances.</li>
<li> Different types of methods to initialize the FiniteVolume.</li>
<li> Different pipe friction models (ConstantLaminar, ConstantTurbulent,
     DetailedFriction) to describe the pressure loss due to the wall friction.</li>
</ul>
</html>"));
  initial equation 
    // Initial conditions
    if initType == NoInit then
      // no initial equations
    elseif initType == InitialValues then
      if not Medium.singleState then
        medium.p = p_start;
      end if;
      if use_T_start then
        medium.T = T_start;
      else
        medium.h = h_start;
      end if;
      medium.Xi = X_start[1:Medium.nXi];
    elseif initType == SteadyState then
      if not Medium.singleState then
         der(medium.p) = 0;
      end if;
      der(medium.h) = 0;
      der(medium.Xi) = zeros(Medium.nXi);
    elseif initType == SteadyStateHydraulic then
      if not Medium.singleState then
         der(medium.p) = 0;
      end if;
      if use_T_start then
        medium.T = T_start;
      else
        medium.h = h_start;
      end if;
      medium.Xi = X_start[1:Medium.nXi];
    else
      assert(false, "Unsupported initialization option");
    end if;
  equation 
    /*
  LongPipes.Components.PipeFriction friction[pipe.n](
    each from_dp=false, 
    each dp_nominal=500/pipe.n, 
    each roughness=1, 
    each diameter=30, 
    each length=length/pipe.n);
*/
      
    // Simple linear pressure drop in each segment
    dp_a = dp_nominal*(if linearPressureDrop then m_flow_a/m_flow_nominal else 
                          abs(m_flow_a)*m_flow_a/m_flow_nominal^2);
    dp_b = dp_nominal*(if linearPressureDrop then -m_flow_b/m_flow_nominal else 
                          abs(-m_flow_b)*(-m_flow_b)/m_flow_nominal^2);
  end PipeSegment;
    
  model PipeFriction 
      "Computes different types of pressure losses in pipes due to friction" 
      
      import FT = Modelica_Fluid.Types.FrictionTypes;
      import CT = Modelica_Fluid.Types.CrossSectionTypes;
      import Modelica.Math;
      
  /* This model requires eta and d as input and provides
   an equation m_flow = f1 (dp) or dp = f2(m_flow)
*/
    input Modelica.SIunits.DynamicViscosity eta 
        "Dummy or upstream dynamic viscosity for detailed friction model used for pressure loss calculation";
    input Modelica.SIunits.Density d 
        "Dummy or upstream density for detailed friction model used for pressure loss calculation";
    Modelica.SIunits.Pressure dp(start=0) "Pressure loss due to pipe friction";
    Modelica.SIunits.MassFlowRate m_flow(start=0) 
        "Mass flow rate from port_a to port_b";
      
    parameter Modelica_Fluid.Types.FrictionTypes.Temp frictionType=Modelica_Fluid.Types.
        FrictionTypes.ConstantTurbulent 
        "Type of friction to determine pressure loss";
    parameter Modelica.SIunits.AbsolutePressure dp_nominal(min=1.e-10)=
        Modelica.SIunits.Conversions.from_bar(1.0) " Nominal pressure drop" 
      annotation (Dialog(enable=frictionType==FT.ConstantLaminar or frictionType==FT.ConstantTurbulent, group=
            "frictionType = ConstantLaminar or ConstantTurbulent"));
      
    parameter Modelica.SIunits.MassFlowRate m_flow_nominal(min=1.e-10)=1 
        " Nominal mass flow rate at nominal pressure drop" 
                                                         annotation (Dialog(
           enable=frictionType==FT.ConstantLaminar or frictionType==FT.ConstantTurbulent, group=
           "frictionType = ConstantLaminar or ConstantTurbulent"));
    parameter Modelica.SIunits.Length length=1 " Length of pipe" 
      annotation (Dialog(enable=frictionType==FT.DetailedFriction, group="frictionType = DetailedFriction"));
    parameter Modelica.SIunits.Length roughness=0 " Roughness of pipe" 
      annotation (Dialog(enable=frictionType==FT.DetailedFriction, group="frictionType = DetailedFriction"));
    parameter Modelica_Fluid.Types.CrossSectionTypes.Temp crossSectionType=
                       Modelica_Fluid.Types.CrossSectionTypes.Circular 
        " Type of cross section of pipe" 
      annotation (Dialog(enable=frictionType==FT.DetailedFriction, group="frictionType = DetailedFriction"));
    parameter Modelica.SIunits.Diameter diameter=0.1 " Inner diameter of pipe" 
      annotation (Dialog(enable=frictionType==FT.DetailedFriction and crossSectionType==CT.Circular, group="frictionType = DetailedFriction"));
    parameter Modelica.SIunits.Length width=0.05 " Inner width of pipe" 
      annotation (Dialog(enable=frictionType==FT.DetailedFriction and crossSectionType==CT.Rectangular, group="frictionType = DetailedFriction"));
    parameter Modelica.SIunits.Length height=0.02 " Inner height of pipe" 
      annotation (Dialog(enable=frictionType==FT.DetailedFriction and crossSectionType==CT.Rectangular, group="frictionType = DetailedFriction"));
    parameter Modelica.SIunits.Area area=0.01 " Cross sectional area of pipe" 
      annotation (Dialog(enable=frictionType==FT.DetailedFriction and crossSectionType==CT.General, group="frictionType = DetailedFriction"));
    parameter Modelica.SIunits.Length perimeter=0.1 
        " Wetted perimeter of cross sectional area" 
      annotation (Dialog(enable=frictionType==FT.DetailedFriction and crossSectionType==CT.General, group="frictionType = DetailedFriction"));
    parameter Boolean from_dp=true 
        " = true, use m_flow = f(dp) otherwise use dp = f(m_flow), i.e., inverse equation"
      annotation (Evaluate=true, Dialog(tab="Advanced"));
    parameter Modelica.SIunits.Pressure p_small(min=1.e-10)=1 
        " A small laminar region is introduced around p_small" 
                                                             annotation (Dialog(
          tab="Advanced", group="Only for frictionType = ConstantTurbulent"));
      
    annotation (
  Images(Parameters(group="frictionType = ConstantLaminar or ConstantTurbulent", source=""),
         Parameters(group="frictionType = DetailedFriction", source="Images/PipeFriction1_small.png")),
  structurallyIncomplete = true,
  preferedView="info",
      Diagram,
      Icon,
      Documentation(info="<html>
<p>
This component models the pressure loss in a short pipe
due to friction under the assumption of quasi steady state flow (i.e., the
mass flow rate varies only slowly). This model is not complete
but may be used in a pipe model to provide an equation to compute
the friction pressure loss from the mass flow rate through
the pipe (see, e.g., <a href=\"Modelica://Modelica_Fluid.Components.ShortPipe\">Modelica_Fluid.Components.ShortPipe</a>).
</p>
<p>
Three loss models can be selected via
parameter <b>frictionType</b>:
</p>
<pre>
   frictionType = <b>ConstantLaminar</b>  :  dp =  k*m_flow
                = <b>ConstantTurbulent</b>:  dp =  k*m_flow^2  if m_flow &gt; 0
                                         = -k*m_flow^2  if m_flow &lt; 0
                = <b>DetailedFriction</b> :  dp = lambda(Re,Delta)*(L*rho/D)*v^2/2
                                         = lambda2(Re,Delta)*L*eta^2/(2*D^3*rho^3)
</pre>
<p>
where dp = \"port_a.p - port_b.p\" is the pressure loss and
m_flow is the mass flow rate from port_a to port_b.
</p>
<h3>ConstantLaminar and ConstantTurbulent</h3>
<p>
The pressure loss factor \"k\" is computed by providing the
mass flow rate \"m_flow_nominal\" and the corresponding
pressure loss \"dp_nominal\" for one flow condition
(usually the desired nominal flow condition). These factors might
be estimated or determined by measurements.
</p>
<p>
For \"ConstantTurbulent\" a small laminar region
is introduced around zero mass flow rate by interpolating
with a cubic polynomial (this technique is copied from the
ThermoFluid library).
</p>
<p>
The first two formulations are useful, if the pipe data is directly
measured and the main operating points are either fully in the
laminar or fully in the turbulent region. It would be better
for \"ConstantTurbulent\" to use the \"real\" laminar region. However,
then more data is required, especially the viscosity and the
diameter of the pipe.
</p>
<h3>DetailedFriction</h3>
<p>
The \"DetailedFriction\" option provides a detailed model
of frictional losses for commercial pipes with
<b>nonuniform roughness</b> (including the smooth pipe
as a special case). For pipes with circular cross section
the pressure loss is computed as:
</p>
<pre>
   dp = lambda*(L/D)*rho*v^2/2
      = lambda2*(L/(2*D^3))*(eta^2/rho)
        (with lambda2 = lambda*Re^2)
</pre>
<p>
where
</p>
<ul>
<li> L is the length of the pipe,</li>
<li> D is the diameter of the pipe,</li>
<li> lambda = lambda(Re,<font face=\"Symbol\">D</font>) is the \"usual\" friction coefficient,</li>
<li> lambda2 = lambda*Re^2 is the friction coefficient used in this model,</li>
<li> Re = v*D*rho/eta is the Reynolds number</li>
<li> <font face=\"Symbol\">D</font> = <font face=\"Symbol\">d</font>/D is the relative roughness where
     \"<font face=\"Symbol\">d</font>\" is
     the absolute \"roughness\", i.e., the averaged height of asperities in the pipe
     (<font face=\"Symbol\">d</font> may change over time due to growth of surface asperities during
      service, see <i>[Idelchick 1994, p. 85, Tables 2-1, 2-2])</i>,</li>
<li> rho is the density,</li>
<li> eta is the dynamic viscosity, and </li>
<li> v is the mean velocity.</li>
</ul>
<p>
The first form is usually given in books but is not suited
for a simulation program since lambda is infinity for zero mass flow rate.
The second form is the one implemented
in this model (lambda2=0 for zero mass flow rate).
The friction coefficient <b>lambda</b> is shown in the next figure:
</p>
<IMG SRC=\"../Images/Components/PipeFriction1.png\" ALT=\"PipeFriction1\">
<p>
More useful for a simulation model is the slightly
differently defined friction coefficient <b>lambda2</b> = lambda*Re^2,
as shown in the next figure:
</p>
<IMG SRC=\"../Images/Components/PipeFriction2.png\" ALT=\"PipeFriction2\">
<p>
<ul>
<li> For <b>Re &le; 2000</b>, the flow is <b>laminar</b> and the exact solution of the
     3-dim. Navier-Stokes equations (momentum and mass balance) is used under the
     assumptions of steady flow, constant pressure gradient and constant
     density and viscosity (= Hagen-Poiseuille flow). </li>
<li> For <b>Re &ge; 4000</b>, the flow is <b>turbulent</b>.
     Depending on the calculation direction (see \"Inverse formulation\"
     below) either of two explicite equations are used. If the pressure drop is assumed
     known (and therefore implicitly also lambda2), then the
     corresponding Reynolds number is computed with the Colebrook-White equation
     <i>[Colebrook 1939; Idelchik 1994, p. 83, eq. (2-9)]</i>.
     These are the <b>red</b> curves in the diagrams above.
     If the mass flow rate is assumed known (and therefore implicitly
     also the Reynolds number), then lambda2 is computed by an approximation of the
     inverse of the Colebrook-White equation <i>[Swamee and Jain 1976;
     Miller 1990, p. 191, eq.(8.4)]</i>.</li>
<li> For <b>2000 &le; Re &le; 4000</b> there is a transition region between laminar
     and turbulent flow. The value of lambda2 depends on more factors as just
     the Reynolds number and the relative roughness, therefore only crude approximations
     are possible in this area.<br>
     The deviation from the laminar region depends on the
     relative roughness. A laminar flow at Re=2000 is only reached for smooth pipes.
     The deviation Reynolds number Re1 is computed according to
     <i>[Samoilenko 1968; Idelchik 1994, p. 81, sect. 2.1.21].</i>
     These are the <b>blue</b> curves in the diagrams above.<br>
     Between Re1=Re1(<font face=\"Symbol\">d</font>/D) and Re2=4000, lambda2 is approximated by a cubic
     polynomial in the \"lg(lambda2) - lg(Re)\" chart (see figure above) such that the
     first derivative is continuous at these two points. In order to avoid
     the solution of non-linear equations, two different cubic polynomials are used
     for the direct and the inverse formulation. This leads to some discrepancies
     in lambda2 (= differences between the red and the blue curves).
     This is acceptable, because the transition region is anyway not
     precisely known since the actual friction coefficient depends on
     additional factors and since the operating points are usually
     not in this region.</li>
</ul>
<p>
The absolute roughness <font face=\"Symbol\">d</font> has usually to
be estimated. In <i>[Idelchik 1994, pp. 105-109,
Table 2-5; Miller 1990, p. 190, Table 8-1]</i> many examples are given.
As a short summary:
</p>
<table border=1 cellspacing=0 cellpadding=2>
  <tr><td><b>Smooth pipes</b></td>
      <td>Drawn brass, coper, aluminium, glass, etc.</td>
      <td><font face=\"Symbol\">d</font> = 0.0025 mm</td>
  </tr>
  <tr><td rowspan=\"3\"><b>Steel pipes</b></td>
      <td>New smooth pipes</td>
      <td><font face=\"Symbol\">d</font> = 0.025 mm</td>
  </tr>
  <tr><td>Mortar lined, average finish</td>
      <td><font face=\"Symbol\">d</font> = 0.1 mm</td>
  </tr>
  <tr><td>Heavy rust</td>
      <td><font face=\"Symbol\">d</font> = 1 mm</td>
  </tr>
  <tr><td rowspan=\"3\"><b>Concrete pipes</b></td>
      <td>Steel forms, first class workmanship</td>
      <td><font face=\"Symbol\">d</font> = 0.025 mm</td>
  </tr>
  <tr><td>Steel forms, average workmanship</td>
      <td><font face=\"Symbol\">d</font> = 0.1 mm</td>
  </tr>
  <tr><td>Block linings</td>
      <td><font face=\"Symbol\">d</font> = 1 mm</td>
  </tr>
</table>
<p>
The equations above are valid for incompressible flow.
They can also be applied for <b>compressible</b> flow up to about <b>Ma = 0.6</b>
(Ma is the Mach number) with a maximum error in lambda of about 3 %.
The effect of gas compressibility in a wide region can be taken into
account by the following formula derived by Voronin
<i>[Voronin 1959; Idelchick 1994, p. 97, sect. 2.1.81]</i>:
</p>
<pre>
  lambda_comp = lambda*(1 + (kappa-1)/2 * Ma^2)^(-0.47)
        kappa = cp/cv // specific heat ratio
</pre>
<p>
An appreciable decrease in the coefficent \"lambda_comp\" is observed
only in a narrow transonic region and also at supersonic flow velocities
by about 15% <i>[Idelchick 1994, p. 97, sect. 2.1.81]</i>.
</p>
<h3>Inverse formulation</h3>
<p>
In the \"Advanced menu\" it is possible via parameter
\"from_dp\" to define in which form the
loss equation is actually evaluated (<b>default</b> is from_dp = <b>true</b>):
</p>
<pre>
   from_dp = <b>true</b>:   m_flow = f1(dp)
           = <b>false</b>:  dp    = f2(m_flow)
</pre>
<p>
\"from_dp\" can be useful to avoid nonlinear systems of equations
in cases where the inverse pressure loss function is needed.
</p>
<p>
At the 34th Modelica meeting in Vienna it was discussed to introduce
a language element for alternatives, such that the tool can
figure out what alternative to use. If this would be available,
parameter from_dp could be removed and the equations would
be written as:
</p>
<pre>
  alternative
    // m_flow = f1(dp);
  or
    // dp = f2(m_flow);
  end alternative;
</pre>
<p>
The tool has then \"somehow\" to select the better alternative.
Further research is needed to develop appropriate symbolic
transformation algorithms.
If you have examples where this is an issue, please provide
them, in order that it is possible to experiment with.
</p>
<h3>References</h3>
<dl><dt>Colebrook F. (1939):</dt>
    <dd><b>Turbulent flow in pipes with particular reference to the transition
         region between the smooth and rough pipe laws</b>.
         J. Inst. Civ. Eng. no. 4, 14-25.</dd>
    <dt>Idelchik I.E. (1994):</dt>
    <dd><a href=\"http://www.begellhouse.com/books/00c0f05b040d2ec0.html\"><b>Handbook
        of Hydraulic Resistance</b></a>. 3rd edition, Begell House, ISBN
        0-8493-9908-4</dd>
    <dt>Miller D. S. (1990):</dt>
    <dd><b>Internal flow systems</b>.
    2nd edition. Cranfield:BHRA(Information Services).</dd>
    <dt>Samoilenko L.A. (1968):</dt>
    <dd><b>Investigation of the Hydraulic Resistance of Pipelines in the
        Zone of Transition from Laminar into Turbulent Motion</b>.
        Thesis (Cand. of Technical Science), Leningrad.</dd>
    <dt>Swamee P.K. and Jain A.K. (1976):</dt>
    <dd><b>Explicit equations for pipe-flow problems</b>.
         Proc. ASCE, J.Hydraul. Div., 102 (HY5), pp. 657-664.</dd>
    <dt>Voronin F.S. (1959):</dt>
    <dd><b>Effect of contraction on the friction coefficient in a
           turbulent gas flow</b>.
           Inzh. Fiz. Zh., vol. 2, no. 11, pp. 81-85.</dd>
</dl>
</html>",   revisions="<html>
<h3>Author</h3>
<p>
<a href=\"http://www.robotic.dlr.de/Martin.Otter/\">Martin Otter</a><br>
Deutsches Zentrum f&uuml;r Luft und Raumfahrt e.V. (DLR)<br>
Institut f&uuml;r Robotik und Mechatronik<br>
Postfach 1116<br>
D-82230 Wessling<br>
Germany<br>
email: <A HREF=\"mailto:Martin.Otter@dlr.de\">Martin.Otter@dlr.de</A><br>
</p>
</html>"));
    Modelica.SIunits.ReynoldsNumber Re 
        "Dummy or Reynolds number of flow, if frictionType = DetailedFriction";
    Real lambda 
        "Dummy or friction coefficient, if frictionType = DetailedFriction";
    Real lambda2 
        "Dummy or non-standard friction coefficient, if frictionType = DetailedFriction (= lambda*Re^2)";
    final parameter Real Delta=roughness/D "Relative roughness";
      
    // Auxiliary variables for ConstantLaminar and ConstantTurbulent
    protected 
    constant Modelica.SIunits.MassFlowRate unitMassFlowRate=1;
    parameter Real k=if frictionType == FT.ConstantLaminar then 
        dp_nominal/m_flow_nominal else (if frictionType == FT.ConstantTurbulent then 
       dp_nominal/m_flow_nominal^2 else length/(2*D*D*D)) 
        "Pressure loss coefficient (dp = k*f(m_flow))";
    parameter Real delta=if from_dp then p_small else sqrt(dp_nominal/k);
    parameter Real C1=if from_dp then 0.5/sqrt(delta) - 3.0*C3*delta^2 else 0.5
        *delta "Coefficient 1 of cubic polynomial in the laminar region";
    parameter Real C3=if from_dp then -0.25/(sqrt(delta)*delta^2) else 0.5/
        delta "Coefficient 3 of cubic polynomial in the laminar region";
      
    // Auxiliary variables for DetailedFriction model
    parameter Modelica.SIunits.Diameter D=if crossSectionType == CT.Circular then 
                diameter else (if crossSectionType == CT.Rectangular then 4*
          width*height/(2*(width + height)) else 4*area/perimeter) 
        "Diameter of pipe in SI units";
    parameter Modelica.SIunits.ReynoldsNumber Re1=(745*exp(if Delta <= 0.0065 then 
                1 else 0.0065/Delta))^(if from_dp then 0.97 else 1) 
        "Re leaving laminar curve";
    parameter Modelica.SIunits.ReynoldsNumber Re2=4000 
        "Re entering turbulent curve";
      
    // point lg(lambda2(Re1)) with derivative at lg(Re1)
    parameter Real x1=if from_dp then Math.log10(64*Re1) else Math.log10(Re1);
    parameter Real y1=if from_dp then Math.log10(Re1) else Math.log10(64*Re1);
    parameter Real yd1=1;
      
    // Point lg(lambda2(Re2)) with derivative at lg(Re2)
    parameter Real aux1=(0.5/Math.log(10))*5.74*0.9;
    parameter Real aux2=Delta/3.7 + 5.74/Re2^0.9;
    parameter Real aux3=Math.log10(aux2);
    parameter Real L2=0.25*(Re2/aux3)^2;
    parameter Real aux4=2.51/sqrt(L2) + 0.27*Delta;
    parameter Real aux5=-2*sqrt(L2)*Math.log10(aux4);
    parameter Real x2=if from_dp then Math.log10(L2) else Math.log10(Re2);
    parameter Real y2=if from_dp then Math.log10(aux5) else Math.log10(L2);
    parameter Real yd2=if from_dp then 0.5 + (2.51/Math.log(10))/(aux5*aux4) else 
              2 + 4*aux1/(aux2*aux3*(Re2)^0.9);
      
    // Constants: Cubic polynomial between lg(Re1) and lg(Re2)
    parameter Real diff_x=x2 - x1;
    parameter Real m=(y2 - y1)/diff_x;
    parameter Real c2=(3*m - 2*yd1 - yd2)/diff_x;
    parameter Real c3=(yd1 + yd2 - 2*m)/(diff_x*diff_x);
    parameter Real lambda2_1=64*Re1;
    constant Real pi=Modelica.Constants.pi;
    Real dx;
    Real aux7;
  equation 
    if frictionType <> FT.DetailedFriction then
      // Assign dummy values for auxiliary variables
      Re = 0;
      dx = 0;
      lambda = 0;
      lambda2 = 0;
      aux7 = 0;
    else
      lambda = noEvent(if Re < 64 then 1 else lambda2/(Re*Re));
    end if;
      
    if from_dp then
      // equations in the form m_flow = m_flow(dp)
      if frictionType == FT.ConstantLaminar then
        m_flow = dp/k;
      elseif frictionType == FT.ConstantTurbulent then
        m_flow = noEvent(if dp > delta then sqrt(dp) else (if dp < -delta then -
          sqrt(-dp) else (C1 + C3*dp*dp)*dp))/sqrt(k);
      else
        lambda2 = noEvent(d*abs(dp)/(k*eta*eta));
        if noEvent(lambda2/64 <= Re1) then
          aux7 = 0;
          dx = 0;
          Re = lambda2/64;
        else
          aux7 = -2*sqrt(lambda2)*Math.log10(2.51/sqrt(lambda2) + 0.27*Delta);
          dx = if noEvent(aux7 >= Re2) then 0 else Math.log10(lambda2/lambda2_1);
          Re = if noEvent(aux7 >= Re2) then aux7 else Re1*(lambda2/lambda2_1)^(
            1 + dx*(c2 + dx*c3));
        end if;
        m_flow = noEvent((pi*D/4)*eta*Re*(if dp >= 0 then +1 else -1));
      end if;
    else
      // equations in the form dp = dp(m_flow)
      if frictionType == FT.ConstantLaminar then
        dp = k*m_flow;
      elseif frictionType == FT.ConstantTurbulent then
        dp = k*noEvent(if m_flow > delta then m_flow*m_flow else (if m_flow < -
          delta then -m_flow*m_flow else (C1 + C3*m_flow*m_flow)*m_flow));
      else
        Re = noEvent((4/pi)*abs(m_flow)/(D*eta));
        dx = noEvent(if Re < Re1 or Re > Re2 then 0 else Math.log10(Re/Re1));
        lambda2 = noEvent(if Re <= Re1 then 64*Re else (if Re >= Re2 then 0.25*
          (Re/Math.log10(Delta/3.7 + 5.74/Re^0.9))^2 else 64*Re1*(Re/Re1)^(1 +
          dx*(c2 + dx*c3))));
        aux7 = 0;
        dp = noEvent(k*lambda2*eta*eta/d*(if m_flow >= 0 then 1 else -1));
      end if;
    end if;
  end PipeFriction;
  end Utilities;
end Modelica_Fluid_WorkInProgress;
