package Volumes "Generic volume, tank and other volume type components" 
   extends Modelica_Fluid.Icons.VariantLibrary;
  
    model MixingVolume 
    "Mixing volume with inlet and outlet ports (flow reversal is allowed)" 
    extends Modelica_Fluid.Interfaces.PartialLumpedVolume(V_lumped=V, Ws_flow=0);
    parameter SI.Volume V "Volume";
    Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a thermalPort 
      "Thermal port" 
      annotation (extent=[-20,88; 20,108]);
    annotation (
      Icon(
        Ellipse(extent=[-100,100; 100,-100], style(
            color=0,
            rgbcolor={0,0,0},
            gradient=3,
            fillColor=68,
            rgbfillColor={170,213,255})),
        Text(extent=[-144,178; 146,116], string="%name"),
        Text(
          extent=[-130,-108; 144,-150],
          style(color=0),
          string="V=%V")),
      Documentation(info="<html>
Ideally mixed volume of constant size with two fluid ports and one medium model. The flow properties are computed from the upstream quantities, pressures are equal in both nodes and the medium model. Heat transfer through a thermal port is possible, it equals zero if the port remains unconnected. The thermal port temperature is equal to the medium temperature.
</html>"),
      Diagram);
    equation 
    thermalPort.T = medium.T;
    Qs_flow = thermalPort.Q_flow;
    end MixingVolume;
  
model OpenTank "Open tank with inlet/outlet ports at the bottom" 
    replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
    annotation (choicesAllMatching=true);
    
//Transformation of kinetic energy
    parameter Boolean p_static_at_port=false 
      "=true, kinetic energy and dissipation is accounted for in port pressure"
                                                                                                        annotation(Evaluate=true, Dialog(tab="Advanced"));
    parameter Real[n_ports] zeta_in=fill(0, n_ports) 
      "Hydraulic resistance into tank, 1 for total dissipation of kinetic energy and uniform flow distribution in pipe"
                                                                                                        annotation(Dialog(tab="Advanced",enable=p_static_at_pot==false));
    parameter Real[n_ports] zeta_out=fill(1, n_ports) 
      "Hydraulic resistance out of tank, 0 for ideal smooth outlet" 
                                                                  annotation(Dialog(tab="Advanced",enable=p_static_at_pot==false));
    
//Tank geometry  
    parameter SI.Height height "Height of tank";
    parameter SI.Area area "Area of tank";
    parameter SI.Volume V0=0 "Volume of the liquid when the level is zero";
    
//Port definitions 
    parameter Integer n_ports(min=1) = 1 "Number of bottom ports (min=1)" 
     annotation(Dialog(group="bottomPorts (= pipes at bottom of tank; in and out flow of tank)"));
    parameter SI.Diameter pipe_diameters[n_ports] 
      "Inner (hydraulic) diameters of bottom ports (array)" 
     annotation(Dialog(group="bottomPorts (= pipes at bottom of tank; in and out flow of tank)", enable=n_bottomPorts > 0));
    Modelica.Fluid.Interfaces.FluidPort_a port[n_ports](
      redeclare package Medium = Medium,
      m_flow(each start=0),
      mXi_flow(each start=0)) 
    annotation (extent=[-12,-109; 8,-89]);
    
//Ambient  
    parameter Medium.AbsolutePressure p_ambient=ambient.default_p_ambient 
      "Tank surface pressure" 
    annotation(Dialog(tab = "Ambient and Initialization", group = "Ambient"));
    parameter Medium.Temperature T_ambient=ambient.default_T_ambient 
      "Tank surface Temperature" 
    annotation(Dialog(tab = "Ambient and Initialization", group = "Ambient"));
    outer Modelica_Fluid.Ambient ambient;
    
//Initialization
    parameter Types.Init.Temp initType=Types.Init.InitialValues 
      "Initialization option" 
    annotation(Evaluate=true,Dialog(tab = "Ambient and Initialization", group = "Initialization"));
    parameter SI.Height level_start "Start value of tank level" 
    annotation(Dialog(tab="Ambient and Initialization", group = "Initialization"));
    parameter Medium.Temperature T_start=T_ambient "Start value of temperature"
    annotation(Dialog(tab = "Ambient and Initialization", group = "Initialization"));
    parameter Medium.MassFraction X_start[Medium.nX]=Medium.X_default 
      "Start value of mass fractions m_i/m" 
    annotation (Dialog(tab="Ambient and Initialization", group = "Initialization", enable=Medium.nXi > 0));
    
//Tank properties  
    Medium.BaseProperties medium(
      preferredMediumStates=true,
      p(start=p_ambient),
      T(start=T_start),
      h(start=h_start),
      Xi(start=X_start[1:Medium.nXi]));
    SI.Height level(stateSelect=StateSelect.default, start=level_start) 
      "Level height of tank";
    SI.Volume V(stateSelect=StateSelect.never) "Actual tank volume";
    SI.Energy U "Internal energy of tank volume";
    SI.Mass m "Mass of fluid in tank";
    SI.Mass mXi[Medium.nXi] "Masses of independent components in the fluid";
    SI.Pressure p_static "bottom tank pressure";
    
  protected 
    parameter Medium.SpecificEnthalpy h_start=Medium.specificEnthalpy_pTX(
        p_ambient,
        T_start,
        X_start);
    parameter SI.Area[n_ports] pipeArea=Modelica.Constants.pi/4*{pipe_diameters[
        i]^2 for i in 1:n_ports};
    
equation 
  //Total quantities
    V = area*level + V0 "Volume of fluid";
    m = V*medium.d "Mass of fluid";
    mXi = m*medium.Xi "Mass of fluid components";
    U = m*medium.u "Internal energy of fluid";
    medium.p = p_ambient;
    
  // Mass balances
    der(m) = sum(port.m_flow);
    for i in 1:Medium.nXi loop
      der(mXi[i]) = sum(port[:].mXi_flow[i]);
    end for;
    
  // Energy balance
    if Medium.singleState then
      der(U) = sum(port.H_flow);
                               //Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)
    else
      der(U) = sum(port.H_flow) - p_ambient*der(V);
    end if;
    assert(level <= height, "Tank is full (level = height = " + String(level) + ")");
    assert(level > 0, "Tank is empty (level = 0), tank model is not designed to allow air flow through ports");
    
//Determine port properties  
    p_static = level*ambient.g*medium.d + p_ambient;
    for i in 1:n_ports loop
      port[i].H_flow = semiLinear(
        port[i].m_flow,
        port[i].h,
        medium.h);
      port[i].mXi_flow = semiLinear(
        port[i].m_flow,
        port[i].Xi,
        medium.Xi);
      if p_static_at_port then
        port[i].p = p_static;
      else
       port[i].p = p_static - smooth(2, noEvent(port[i].m_flow^2/(2*medium.d*
          pipeArea[i]^2)*(if port[i].m_flow < 0 then (1 + zeta_out[i]) else (1
           - zeta_in[i]))));
       end if;
    end for;
    
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
    annotation (
      Icon(
        Rectangle(extent=[-100,100; 100,0], style(color=7, fillColor=7)),
        Rectangle(extent=DynamicSelect([-100,-100; 100,10], [-100,-100; 100,(-100
               + 200*level/height)]), style(
            color=69,
            rgbcolor={0,127,255},
            fillColor=71,
            rgbfillColor={85,170,255},
            fillPattern=1)),
        Line(points=[-100,100; -100,-100; 100,-100; 100,100], style(
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
          extent=[-129,53; 130,39],
          style(color=0),
          Line(points=[-100,100; 100,100], style(
              color=0,
              rgbcolor={0,0,0},
              pattern=3)),
        string="start = %level_start m"),
        Text(
          extent=[-95,30; 95,5],
          style(color=0),
          string=DynamicSelect(" ", realString(
              level,
              1,
              integer(precision)))),
        Line(points=[-100,100; 100,100], style(
            color=0,
            rgbcolor={0,0,0},
            pattern=3))),
      Documentation(info="<HTML>
<p>
This is a simplified model of a tank. 
The top part is open to the environment at the fixed pressure 
<tt>p_ambient</tt>. Heat transfer to the environment and to 
the tank walls is neglected.
The tank is filled with a single or multiple-substance liquid, 
assumed to have uniform temperature and mass fractions.</p> 
<p>
Inlet and outlet connections are situated at the bottom of the tank. The following assumptions are made:
</p>
<ul>
<li>Incompressible medium, liquid fluid is of uniform density and uniform temperature</li>
<li>Heat transfer to the environment is neglected</li>
<li>Kinetic energy of the fluid in the tank is neglected, the cross sectional area of the tank is larger than the cross sectional area of the inlet/outlet by several orders of magnitude</li>
<li>No air is leaving the tank through the ports, if the liquid level drops below zero the simulation stops.</li>
</ul>
<p>By default the port pressure is the pressure just after the outlet (or just before the inlet) in the attached pipe. The hydraulic resistances <tt>zeta_in</tt> and <tt>zeta_out</tt> determine the dissipative pressure drop between tank and port depending on the direction of mass flow. The default values (zeta_in=1, zeta_out=0) assume no dissipation at the tank outlet (ideal smooth opening) and total dissipation of kinetic energy at the tank inlet. Larger values are found for sharp edged openings and non-uniform velocity distributions in the pipe. A large selection of possible cases are listed in <i>[Idelchik, Handbook of Hydraulic Resistance, 2004]</i>. If the flag <tt>static_pressure_at_port</tt> in the <tt>Advanced</tt> menu is set to true, the port pressure represents the static head at the bottom of the tank. The relationship between pressure drop and mass flow rate must then be provided by the connected component. </p>  
 
 
</HTML>", revisions="<html>
<ul>
<li><i>Jan. 6, 2006</i> by Katja Poschlad, Manuel Remelhe (AST Uni Dortmund), 
   Martin Otter (DLR):<br> 
   Implementation based on former tank model.</li>
<li><i>Apr. 25, 2006</i> by Katrin Pr&ouml;l&szlig; (TUHH):<br>
Limitation to bottom ports only, added inlet and outlet loss factors.</li>
</ul>
</html>"),
      Diagram,
      uses(Modelica(version="2.2.1"), Modelica_Fluid(version="0.952")),
      Coordsys(grid=[1,1], scale=0.2));
end OpenTank;
  
model Tank 
    "Open tank with top and bottom inlet/outlet ports at a defineable height" 
    
  import Modelica.Constants;
  import Modelica_Fluid.PressureLosses.BaseClasses.lossConstant_D_zeta;
  import Modelica_Fluid.Utilities.regRoot2;
  import Modelica_Fluid.Volumes.BaseClasses.TankPortData;
    
  replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
    annotation (choicesAllMatching=true);
    
  SI.Height level(stateSelect=StateSelect.prefer, start=level_start) 
      "Fluid level in the tank";
    
//Tank geometry  
    parameter SI.Height levelMax "Maximum level of tank before it overflows";
    parameter SI.Area area "Area of tank";
    parameter SI.Volume V0=0 "Volume of the liquid when level = 0";
    
//Port definitions 
    parameter Integer nTopPorts(min=1) = 1 
      "Number of inlet ports above levelMax (>= 1)";
    
    Modelica.Fluid.Interfaces.FluidPorts_a topPorts[nTopPorts](
    redeclare package Medium = Medium,
    m_flow(each start=0, each min=0),
    mXi_flow(each start=0, each min=0)) 
      "Inlet ports over levelMax at top of tank (fluid flows only from the port in to the tank)"
    annotation (extent=[-5,85; 5,125],  rotation=90);
    
    parameter Modelica_Fluid.Volumes.BaseClasses.TankPortData portsData[:] = {TankPortData(diameter=0)} 
      "Data of inlet/outlet ports at side and bottom of tank";
    
    Modelica.Fluid.Interfaces.FluidPorts_b ports[size(portsData,1)](
    redeclare package Medium = Medium,
    m_flow(each start=0),
    mXi_flow(each start=0)) 
      "inlet/outlet ports at bottom or side of tank (fluid flows in to or out of port; a port might be above the fluid level)"
    annotation (extent=[-5,-125; 5,-85],  rotation=90);
    
//Ambient  
   outer Modelica_Fluid.Ambient ambient "Ambient conditions";
   parameter Medium.AbsolutePressure p_ambient=ambient.default_p_ambient 
      "Tank surface pressure" 
    annotation(Dialog(tab = "Ambient and Initialization", group = "Ambient"));
   parameter Medium.Temperature T_ambient=ambient.default_T_ambient 
      "Tank surface Temperature" 
    annotation(Dialog(tab = "Ambient and Initialization", group = "Ambient"));
    
//Initialization
    parameter Types.Init.Temp initType=Types.Init.InitialValues 
      "Initialization option" 
    annotation(Evaluate=true,Dialog(tab = "Ambient and Initialization", group = "Initialization"));
    parameter SI.Height level_start(min=0) "Start value of tank level" 
    annotation(Dialog(tab="Ambient and Initialization", group = "Initialization"));
    parameter Medium.Temperature T_start=T_ambient "Start value of temperature"
    annotation(Dialog(tab = "Ambient and Initialization", group = "Initialization"));
    parameter Medium.MassFraction X_start[Medium.nX]=Medium.X_default 
      "Start value of mass fractions m_i/m" 
    annotation (Dialog(tab="Ambient and Initialization", group = "Initialization", enable=Medium.nXi > 0));
    
// Advanced  
    parameter Real hysteresisFactor(min=0) = 0.1 
      "Hysteresis for empty pipe = diameter*hysteresisFactor" 
    annotation(Dialog(tab="Advanced", group="Numerical properties"));
    parameter SI.MassFlowRate m_flow_small(min=0) = 1e-5 
      "Regularization range at zero mass flow rate" 
    annotation(Dialog(tab="Advanced", group="Numerical properties"));
    parameter Boolean stiffCharacteristicForEmptyPort = true 
      "=true, if steep pressure loss characteristic for empty pipe port" 
    annotation(Dialog(tab="Advanced", group="Numerical properties"), Evaluate=true);
    parameter Real zetaLarge(min=0) = 1e5 
      "Large pressure loss factor if mass flows out of empty pipe port" 
    annotation(Dialog(tab="Advanced", group="Numerical properties", enable=stiffCharacteristicForEmptyPort));
    
//Tank properties  
     final parameter Integer nPorts = size(ports,1) 
      "Number of inlet/outlet ports";
     final parameter Medium.SpecificEnthalpy h_start=Medium.specificEnthalpy_pTX(
        p_ambient,
        T_start,
        X_start) annotation(Hide=true);
    Medium.BaseProperties medium(
      preferredMediumStates=true,
      p(start=p_ambient),
      T(start=T_start),
      h(start=h_start),
      Xi(start=X_start[1:Medium.nXi]));
    SI.Volume V(stateSelect=StateSelect.never) "Actual tank volume";
    SI.Energy U(stateSelect=StateSelect.never) "Internal energy of tank volume";
    SI.Mass m(stateSelect=StateSelect.never) "Mass of fluid in tank";
    SI.Mass mXi[Medium.nXi](each stateSelect=StateSelect.never) 
      "Masses of independent components in the fluid";
    
  protected 
    parameter SI.Area bottomArea[nPorts]=Constants.pi*{(portsData[i].diameter/2)^2 for i in 1:nPorts};
    parameter SI.Diameter ports_emptyPipeHysteresis[nPorts] = portsData.diameter*hysteresisFactor;
    SI.Length levelAbovePort[nPorts] "Height of fluid over bottom ports";
    Boolean ports_m_flow_out[nPorts](each start = true, each fixed=true);
    Boolean aboveLevel[nPorts] "= true, if level >= ports[i].portLevel";
    Real zeta_out[nPorts];
equation 
  assert(level <= levelMax, "Tank starts to overflow (level = levelMax = " + String(level) + ")");
  assert(m>=0, "Mass in tank is zero");
    
  // Total quantities
    medium.p = p_ambient;
    V = area*level + V0 "Volume of fluid";
    m = V*medium.d "Mass of fluid";
    mXi = m*medium.Xi "Mass of fluid components";
    U = m*medium.u "Internal energy of fluid";
    
  // Mass balances
    der(m) = sum(topPorts.m_flow) + sum(ports.m_flow);
    for i in 1:Medium.nXi loop
      der(mXi[i]) = sum(topPorts.mXi_flow[i]) + sum(ports.mXi_flow[i]);
    end for;
    
  // Energy balance
    if Medium.singleState then
      der(U) = sum(topPorts.H_flow) + sum(ports.H_flow);
                               //Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)
    else
      der(U) = sum(topPorts.H_flow) + sum(ports.H_flow) - p_ambient*der(V);
    end if;
    
  // Properties at top ports
    for i in 1:nTopPorts loop
       // It is assumed that fluid flows only into one of the top ports and never out of it 
       topPorts[i].H_flow   = semiLinear(topPorts[i].m_flow, topPorts[i].h, h_start);
       topPorts[i].mXi_flow = semiLinear(topPorts[i].m_flow, topPorts[i].Xi, X_start[1:Medium.nXi]);
       topPorts[i].p        = p_ambient;
/*
       assert(topPorts[i].m_flow > -1, "Mass flows out of tank via topPorts[" + String(i) + "]\n" +
                                         "This indicates a wrong model");
*/
    end for;
    
  // Properties at bottom ports
    for i in 1:nPorts loop
       ports[i].H_flow = semiLinear(ports[i].m_flow, ports[i].h, medium.h);
       ports[i].mXi_flow = semiLinear(ports[i].m_flow, ports[i].Xi, medium.Xi);
       aboveLevel[i] = level >= (portsData[i].portLevel + ports_emptyPipeHysteresis[i])
                       or pre(aboveLevel[i]) and level >= (portsData[i].portLevel - ports_emptyPipeHysteresis[i]);
       levelAbovePort[i] = if aboveLevel[i] then level - portsData[i].portLevel else 0;
      
       if stiffCharacteristicForEmptyPort then
          // If port is above fluid level, use large zeta if fluid flows out of port (= small mass flow rate)
          zeta_out[i] = 1 + (if aboveLevel[i] then 0 else zetaLarge);
          ports[i].p = p_ambient + levelAbovePort[i]*ambient.g*medium.d
                               + Modelica_Fluid.Utilities.regSquare2(ports[i].m_flow, m_flow_small,
                                     lossConstant_D_zeta(portsData[i].diameter, 0.01)/medium.d,
                                     lossConstant_D_zeta(portsData[i].diameter, zeta_out[i])/medium.d);
          ports_m_flow_out[i] = false;
        
       else
          // Handling according to Remelhe/Poschlad
          ports_m_flow_out[i] = (pre(ports_m_flow_out[i]) and not ports[i].p>p_ambient)
                                     or ports[i].m_flow < -1e-6;
         if aboveLevel[i] then
             ports[i].p = p_ambient + levelAbovePort[i]*ambient.g*medium.d -
                               smooth(2,noEvent(if ports[i].m_flow < 0 then ports[i].m_flow^2/
                                     (2*medium.d*bottomArea[i]^2) else 0));
         else
            if pre(ports_m_flow_out[i]) then
               ports[i].m_flow = 0;
            else
               ports[i].p = p_ambient;
            end if;
         end if;
          zeta_out[i] =0;
       end if;
     end for;
    
initial equation 
    for i in 1:nPorts loop
       pre(aboveLevel[i]) = level_start >= portsData[i].portLevel;
    end for;
    
    if initType == Types.Init.NoInit then
    // no initial equations
    elseif initType == Types.Init.InitialValues then
      level = level_start;
      medium.T = T_start;
      medium.Xi = X_start[1:Medium.nXi];
    elseif initType == Types.Init.SteadyState then
      der(level) = 0;
      der(medium.T) = 0;
      der(medium.Xi) = zeros(Medium.nXi);
    elseif initType == Types.Init.SteadyStateHydraulic then
      der(level) = 0;
      medium.T = T_start;
      medium.Xi = X_start[1:Medium.nXi];
    else
      assert(false, "Unsupported initialization option");
    end if;
    
    annotation (
      Icon(
        Rectangle(extent=[-100,-100; 100,100], style(
            color=7,
            rgbcolor={255,255,255},
            fillColor=7,
            rgbfillColor={255,255,255})),
          Rectangle(extent=DynamicSelect([-100,-100; 100,0], [-100,-100; 100,(-100
                 + 200*level/levelMax)]), style(
              color=69,
              rgbcolor={0,127,255},
              fillColor=71,
              rgbfillColor={85,170,255},
              fillPattern=1)),
        Text(
          extent=[-94,19; 96,-1],
        string=DynamicSelect(" ", realString(level, 1, 3)),
        style(color=0, rgbcolor={0,0,0})),
        Line(points=[-100,100; 100,100], style(
            color=0,
            rgbcolor={0,0,0},
            pattern=3)),
      Text(
          extent=[-94,90;95,60],
          style(color=3, rgbcolor={0,0,255}),
          string="%name"),
      Text(
        extent=[-95,-85; 95,-65],
        style(color=0),
          string="%level_start"),
      Text(
        extent=[-95,-55; 95,-35],
        style(color=0),
          string="level_start ="),
      Text(extent=[-95,50; 95,30], string="level =",
        style(color=0, rgbcolor={0,0,0})),
        Line(points=[-100,100; -100,-100; 100,-100; 100,100], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=69,
            rgbfillColor={0,127,255},
            fillPattern=1))),
      Documentation(info="<HTML>
<p> 
Model of a tank that is open to the environment at the fixed pressure
<tt>p_ambient</tt>. Heat transfer to the environment and to 
the tank walls is neglected.
The tank is filled with a single or multiple-substance liquid, 
assumed to have uniform temperature and mass fractions.
</p> 
 
<p>
At the top of the tank over the maximal fill level <b>levelMax</b> 
a vector of FluidPorts, called <b>topPorts</b>, is present.
The assumption is made that fluid flows always in to the tank via these
ports (and never back in to the connector).
If the tank has no top ports, set <b>nTopPorts</b> = 1, and do not
connect to this port (the default connection semantics of Modelica
leads to a behaviour as if the port would not be present; 
the reason is that some tools do currently no support zero sized
connectors).
</p>
 
<p>
The vector of connectors <b>ports</b> are fluid ports at the bottom
and side of the tank at a defineable height. Fluid can flow either out
of or in to this port. The fluid level of the tank may be below
one of these ports. This case is approximated by introducing a
large pressure flow coefficient so that the mass flow rate
through this port is very small in this case.
</p>
 
<p>
If the tank starts to over flow (i.e., level > levelMax), an
assertion is triggered.
</p>
 
<p>
When the diagram layer is open in the plot environment, the
level of the tank is dynamically visualized. Note, the speed
of the diagram animation in Dymola can be set via command
<b>animationSpeed</b>(), e.g., animationSpeed(speed = 10)
</p>
</HTML>", revisions="<html>
<ul>
<li><i>Jul. 29, 2006</i> by Martin Otter (DLR):<br> 
   Improved handling of ports that are above the fluid level and
   simpler implementation.</li>
 
<li><i>Jan. 6, 2006</i> by Katja Poschlad, Manuel Remelhe (AST Uni Dortmund), 
   Martin Otter (DLR):<br> 
   Implementation based on former tank model but with several improvements
   (top, bottom, side ports; correctly treating kinetic energy for outlet
   and total dissipation for inlet; ports can be above the fluid level).</li>
</ul>
</html>"),
      Diagram,
      uses(Modelica(version="2.2.1"), Modelica_Fluid(version="0.952")),
      Coordsys(grid=[1,1], scale=0.2));
end Tank;
  
  package BaseClasses 
    extends Modelica_Fluid.Icons.BaseClassLibrary;
    
    record TankPortData 
      "Data to describe inlet/outlet pipes at the bottom or side of the tank" 
          extends Modelica.Icons.Record;
      
      parameter SI.Diameter diameter 
        "Inner (hydraulic) diameter of inlet/outlet port";
      parameter SI.Height portLevel=0 
        "level of inlet/outlet port (height over the tank base)";
    end TankPortData;
    
  end BaseClasses;
  annotation (Documentation(info="<html>
 
</html>"));
end Volumes;
