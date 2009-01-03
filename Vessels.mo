within Modelica_Fluid;
package Vessels "Devices for storing fluid"
   extends Modelica_Fluid.Icons.VariantLibrary;

    model Volume "Fixed volume with ports, closed to the environment"
      import Modelica.Constants.pi;
      extends Modelica_Fluid.Vessels.BaseClasses.PartialLumpedVolumePorts(
        heatTransfer(surfaceAreas={4*pi*(3/4*V/pi)^(2/3)}));

      parameter SI.Volume V "Volume";

    equation
      fluidVolume = V;
      Ws_flow = 0;
      ports_p_static = medium.p;

      annotation (defaultComponentName="volume",
        Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{
              100,100}}), graphics={Ellipse(
            extent={{-100,100},{100,-100}},
            lineColor={0,0,0},
            fillPattern=FillPattern.Sphere,
            fillColor={170,213,255}), Text(
            extent={{-150,12},{150,-18}},
            lineColor={0,0,0},
            textString="V=%V")}),
      Documentation(info="<html>
Ideally mixed volume of constant size with two fluid ports and one medium model. 
The flow properties are computed from the upstream quantities, pressures are equal in both nodes and the medium model. 
Heat transfer through a thermal port is possible, it equals zero if the port remains unconnected. 
A spherical shape is assumed for the heat transfer area, with V=4/3*pi*r^3, A=4*pi*r^2.
Ideal heat transfer is assumed per default; the thermal port temperature is equal to the medium temperature.
</html>"),
      Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},{
              100,100}}),
              graphics));
    end Volume;

model OpenTank "Open tank with inlet/outlet ports at the bottom"
  import Modelica.Constants.pi;
  extends BaseClasses.PartialLumpedVolumePorts(
    heatTransfer(surfaceAreas={crossArea+2*sqrt(crossArea*pi)*level}),
    final initialize_p = false,
    final p_start = p_ambient,
    final use_d_nominal = false,
    final d_nominal = 0);

  // Tank geometry
  parameter SI.Height height "Height of tank";
  parameter SI.Area crossArea "Area of tank";
  parameter SI.Volume V0=0 "Volume of the liquid when the level is zero";

  // Ambient
  parameter Medium.AbsolutePressure p_ambient=system.p_ambient
      "Tank surface pressure" 
    annotation(Dialog(tab = "Advanced", group = "Ambient"));
  parameter Medium.Temperature T_ambient=system.T_ambient
      "Tank surface Temperature" 
    annotation(Dialog(tab = "Advanced", group = "Ambient"));

  // Initialization
  parameter SI.Height level_start "Start value of tank level" 
    annotation(Dialog(tab="Initialization"));

  // Tank properties
  SI.Height level(stateSelect=StateSelect.default, start=level_start)
      "Level height of tank";
  SI.Volume V(stateSelect=StateSelect.never) "Actual tank volume";

equation
  // Total quantities
  V = crossArea*level + V0 "Volume of fluid";
  fluidVolume = V;
  medium.p = p_ambient;

  // Source termsEnergy balance
  if Medium.singleState or energyDynamics == Types.Dynamics.SteadyState then
    Ws_flow = 0
        "Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)";
  else
    Ws_flow = -p_ambient*der(V);
  end if;

  assert(level <= height, "Tank is full (level = height = " + String(level) + ")");
  assert(level > 0, "Tank is empty (level = 0), tank model is not designed to allow air flow through ports");

  //Determine port properties
  ports_p_static = level*system.g*medium.d + p_ambient;

initial equation
    if massDynamics == Types.Dynamics.FixedInitial then
      level = level_start;
    elseif massDynamics == Types.Dynamics.SteadyStateInitial then
      der(level) = 0;
    end if;

    annotation (
      Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1},
          initialScale=0.2), graphics={
          Rectangle(
            extent={{-100,100},{100,0}},
            lineColor={255,255,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent=DynamicSelect({{-100,-100},{100,10}}, {{-100,-100},{100,(-100
                 + 200*level/height)}}),
            lineColor={0,127,255},
            fillColor={85,170,255},
            fillPattern=FillPattern.Solid),
          Line(points={{-100,100},{-100,-100},{100,-100},{100,100}}, color={0,0,
                0}),
          Text(
            extent={{-129,40},{130,26}},
            lineColor={0,0,0},
            textString="%level_start m"),
          Text(
            extent={{-95,30},{95,5}},
            lineColor={0,0,0},
            textString=DynamicSelect(" ", realString(
                level,
                1,
                integer(precision)))),
          Line(
            points={{-100,100},{100,100}},
            color={0,0,0},
            pattern=LinePattern.Dot),
          Text(
            extent={{-126,81},{133,67}},
            lineColor={0,0,0},
            textString="level_start =")}),
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
<li>No liquid is leaving the tank through the open top, if the liquid level growth over the height the simulation stops.</li>
</ul>
<p>By default the port pressure is the pressure just after the outlet (or just before the inlet) in the attached pipe. The hydraulic resistances <tt>zeta_in</tt> and <tt>zeta_out</tt> determine the dissipative pressure drop between tank and port depending on the direction of mass flow. The default values (zeta_in=1, zeta_out=0) assume no dissipation at the tank outlet (ideal smooth opening) and total dissipation of kinetic energy at the tank inlet. Larger values are found for sharp edged openings and non-uniform velocity distributions in the pipe. A large selection of possible cases are listed in <i>[Idelchik, Handbook of Hydraulic Resistance, 2004]</i>. If the flag <tt>static_pressure_at_port</tt> in the <tt>Advanced</tt> menu is set to true, the port pressure represents the static head at the bottom of the tank. The relationship between pressure drop and mass flow rate must then be provided by the connected component. </p>  
 
 
</HTML>", revisions="<html>
<ul>
<li><i>Dec. 12, 2008</i> by Ruediger Franke: replace energy and mass balance as well as port definitions 
   with common definition in BaseClasses.PartialLumpedVolumePorts</li>
<li><i>Jan. 6, 2006</i> by Katja Poschlad, Manuel Remelhe (AST Uni Dortmund), 
   Martin Otter (DLR):<br> 
   Implementation based on former tank model.</li>
<li><i>Apr. 25, 2006</i> by Katrin Pr&ouml;l&szlig; (TUHH):<br>
Limitation to bottom ports only, added inlet and outlet loss factors.</li>
<li><i>Oct. 29, 2007</i> by Carsten Heinrich (ILK Dresden):<br>
Adapted to the new fluid library interfaces: 
<ul> <li>FluidPorts_b is used instead of FluidPort_b (due to it is defined as an array of ports)</li>
    <li>Port name changed from port to ports</li></ul>Updated documentation.</li>
<li><i>Dec. 8, 2008</i> by Michael Wetter (LBNL):<br>
Implemented trace substances.</li>
</ul>
</html>"),
      Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1},
          initialScale=0.2), graphics),
      uses(Modelica(version="2.2.1"), Modelica_Fluid(version="0.952")));
equation

end OpenTank;

model Tank
    "Open tank with top and bottom inlet/outlet ports at a defineable height"
  extends BaseClasses.PartialLumpedVolume(
    final initialize_p = false,
    final p_start = p_ambient,
    final use_d_nominal = false,
    final d_nominal = 0);

    import Modelica.Constants;
    import Modelica_Fluid.Fittings.BaseClasses.lossConstant_D_zeta;
    import Modelica_Fluid.Utilities.regRoot2;
    import Modelica_Fluid.Vessels.BaseClasses.TankPortData;

  SI.Height level(stateSelect=StateSelect.prefer, start=level_start)
      "Fluid level in the tank";

  //Tank geometry
  parameter SI.Height levelMax "Maximum level of tank before it overflows";
  parameter SI.Area crossArea "Area of tank";
  parameter SI.Volume V0=0 "Volume of the liquid when level = 0";

  //Port definitions
  parameter Integer nTopPorts(min=1) = 1
      "Number of inlet ports above levelMax (>= 1)";

  Modelica_Fluid.Interfaces.FluidPorts_a topPorts[nTopPorts](
    redeclare package Medium = Medium,
    m_flow(each start=0, each min=0))
      "Inlet ports over levelMax at top of tank (fluid flows only from the port in to the tank)"
    annotation (Placement(transformation(
        extent={{0,-20},{10,20}},
        rotation=90,
        origin={0,100})));

  parameter Modelica_Fluid.Vessels.BaseClasses.TankPortData portsData[:]={
        TankPortData(diameter=0.0001)}
      "Data of inlet/outlet ports at side and bottom of tank";

  Modelica_Fluid.Interfaces.FluidPorts_b ports[size(portsData,1)](
    redeclare package Medium = Medium,
    m_flow(each start=0))
      "inlet/outlet ports at bottom or side of tank (fluid flows in to or out of port; a port might be above the fluid level)"
    annotation (Placement(transformation(
        extent={{0,-20},{-10,20}},
        rotation=90,
        origin={0,-100})));

  //Initialization
  parameter SI.Height level_start(min=0) "Start value of tank level" 
    annotation(Dialog(tab="Initialization"));

  // Heat transfer through boundary
  parameter Boolean use_HeatTransfer = false
      "= true to use the HeatTransfer model" 
      annotation (Dialog(tab="Assumptions", group="Heat transfer"));
  replaceable model HeatTransfer = 
      BaseClasses.HeatTransfer.IdealHeatTransfer 
    constrainedby BaseClasses.HeatTransfer.PartialVesselHeatTransfer
      "Wall heat transfer" 
      annotation (Dialog(tab="Assumptions", group="Heat transfer",enable=use_HeatTransfer),choicesAllMatching=true);
  HeatTransfer heatTransfer(
    redeclare final package Medium = Medium,
    n=1,
    states = {medium.state},
    surfaceAreas={crossArea+2*sqrt(crossArea*Modelica.Constants.pi)*level}) 
      annotation (Placement(transformation(
        extent={{-10,-10},{30,30}},
        rotation=90,
        origin={-50,-10})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort if use_HeatTransfer 
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

  // Advanced
  parameter Real hysteresisFactor(min=0) = 0.1
      "Hysteresis for empty pipe = diameter*hysteresisFactor" 
    annotation(Dialog(tab="Advanced", group="Port properties"));
  parameter Boolean stiffCharacteristicForEmptyPort = false
      "=true, if steep pressure loss characteristic for empty pipe port" 
    annotation(Dialog(tab="Advanced", group="Port properties"), Evaluate=true);
  parameter Real zetaLarge(min=0) = 1e5
      "Large pressure loss factor if mass flows out of empty pipe port" 
    annotation(Dialog(tab="Advanced", group="Port properties", enable=stiffCharacteristicForEmptyPort));
  parameter SI.MassFlowRate m_flow_small(min=0) = 1e-5
      "Regularization range at zero mass flow rate" 
    annotation(Dialog(tab="Advanced", group="Port properties", enable=stiffCharacteristicForEmptyPort));

  //Ambient
  parameter Medium.AbsolutePressure p_ambient=system.p_ambient
      "Tank surface pressure" 
    annotation(Dialog(tab = "Advanced", group = "Ambient"));
  parameter Medium.Temperature T_ambient=system.T_ambient
      "Tank surface Temperature" 
    annotation(Dialog(tab = "Advanced", group = "Ambient"));

  // Tank properties
  final parameter Integer nPorts = size(ports,1) "Number of inlet/outlet ports";
  SI.Volume V(stateSelect=StateSelect.never) "Actual tank volume";
  Medium.EnthalpyFlowRate H_flow_top[nTopPorts]
      "Enthalpy flow rates from the top ports in to the tank";
  Medium.EnthalpyFlowRate port_b_H_flow_bottom[nPorts]
      "Enthalpy flow rates from the bottom ports in to the tank";
  Medium.MassFlowRate mXi_flow_top[nTopPorts, Medium.nXi]
      "Substance mass flow rates from the top ports into the tank";
  Medium.MassFlowRate port_b_mXi_flow_bottom[nPorts, Medium.nXi]
      "Substance mass flow rates from the bottom ports into the tank";
  Medium.MassFlowRate mC_flow_top[nTopPorts, Medium.nC]
      "Trace substance mass flow rates from the top ports into the tank";
  Medium.MassFlowRate port_b_mC_flow_bottom[nPorts, Medium.nC]
      "Trace substance mass flow rates from the bottom ports into the tank";
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

  // Only one connection allowed to a port to avoid unwanted ideal mixing
/*
for i in 1:nTopPorts loop
  assert(cardinality(topPorts[i]) <= 1,"
topPorts[" + String(i) + "] of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller.
");
end for;
 
for i in 1:nPorts loop
  assert(cardinality(ports[i]) <= 1,"
ports[" + String(i) + "] of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections which is usually not the intention
of the modeller.
");
end for;
*/

  // Total quantities
  medium.p = p_ambient;
  V = crossArea*level + V0 "Volume of fluid";
  fluidVolume = V;

  // Mass balances
  ms_flow = sum(topPorts.m_flow) + sum(ports.m_flow);
  for i in 1:Medium.nXi loop
    msXi_flow[i] = sum(mXi_flow_top[:,i]) + sum(port_b_mXi_flow_bottom[:,i]);
  end for;
  for i in 1:Medium.nC loop
    msC_flow[i]  = sum(mC_flow_top[:,i])  + sum(port_b_mC_flow_bottom[:,i]);
  end for;

  // Energy balance
  Hs_flow = sum(H_flow_top) + sum(port_b_H_flow_bottom);
  Qs_flow = heatTransfer.Q_flows[1];
  if Medium.singleState or energyDynamics == Types.Dynamics.SteadyState then
    Ws_flow = 0
        "Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)";
  else
    Ws_flow = -p_ambient*der(V);
  end if;

  // Properties at top ports
    for i in 1:nTopPorts loop
       // It is assumed that fluid flows only from one of the top ports in to the tank and never vice versa
       H_flow_top[i]     = topPorts[i].m_flow*actualStream(topPorts[i].h_outflow);
       mXi_flow_top[i,:] = topPorts[i].m_flow*actualStream(topPorts[i].Xi_outflow);
       mC_flow_top[i,:]  = topPorts[i].m_flow*actualStream(topPorts[i].C_outflow);
       topPorts[i].p     = p_ambient;
       topPorts[i].h_outflow = h_start;
       topPorts[i].Xi_outflow = X_start[1:Medium.nXi];
       topPorts[i].C_outflow  = C_start;
/*
       assert(topPorts[i].m_flow > -1, "Mass flows out of tank via topPorts[" + String(i) + "]\n" +
                                         "This indicates a wrong model");
*/
    end for;

  // Properties at bottom ports
    for i in 1:nPorts loop
       port_b_H_flow_bottom[i]   = ports[i].m_flow*actualStream(ports[i].h_outflow);
       port_b_mXi_flow_bottom[i,:] = ports[i].m_flow*actualStream(ports[i].Xi_outflow);
       port_b_mC_flow_bottom[i,:]  = ports[i].m_flow*actualStream(ports[i].C_outflow);
       aboveLevel[i] = level >= (portsData[i].portLevel + ports_emptyPipeHysteresis[i])
                       or pre(aboveLevel[i]) and level >= (portsData[i].portLevel - ports_emptyPipeHysteresis[i]);
       levelAbovePort[i] = if aboveLevel[i] then level - portsData[i].portLevel else 0;
       ports[i].h_outflow = medium.h;
       ports[i].Xi_outflow = medium.Xi;
       ports[i].C_outflow  = C;

       if stiffCharacteristicForEmptyPort then
          // If port is above fluid level, use large zeta if fluid flows out of port (= small mass flow rate)
          zeta_out[i] = 1 + (if aboveLevel[i] then 0 else zetaLarge);
          ports[i].p = p_ambient + levelAbovePort[i]*system.g*medium.d
                               + Modelica_Fluid.Utilities.regSquare2(ports[i].m_flow, m_flow_small,
                                     lossConstant_D_zeta(portsData[i].diameter, 0.01)/medium.d,
                                     lossConstant_D_zeta(portsData[i].diameter, zeta_out[i])/medium.d);
          ports_m_flow_out[i] = false;

       else
          // Handling according to Remelhe/Poschlad
          ports_m_flow_out[i] = (pre(ports_m_flow_out[i]) and not ports[i].p>p_ambient)
                                     or ports[i].m_flow < -1e-6;
         if aboveLevel[i] then
             ports[i].p = p_ambient + levelAbovePort[i]*system.g*medium.d -
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

    if massDynamics == Types.Dynamics.FixedInitial then
      level = level_start;
    elseif massDynamics == Types.Dynamics.SteadyStateInitial then
      der(level) = 0;
    end if;

    annotation (
      Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1},
          initialScale=0.2), graphics={
          Rectangle(
            extent={{-100,-100},{100,100}},
            lineColor={255,255,255},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Rectangle(
            extent=DynamicSelect({{-100,-100},{100,0}}, {{-100,-100},{100,(-100
                 + 200*level/levelMax)}}),
            lineColor={0,127,255},
            fillColor={85,170,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-94,19},{96,-1}},
            lineColor={0,0,0},
            textString=DynamicSelect(" ", realString(
                level,
                1,
                3))),
          Line(
            points={{-100,100},{100,100}},
            color={0,0,0},
            pattern=LinePattern.Dot),
          Text(
            extent={{-94,90},{95,60}},
            lineColor={0,0,255},
            textString="%name"),
          Text(
            extent={{-95,-85},{95,-65}},
            lineColor={0,0,0},
            textString="%level_start"),
          Text(
            extent={{-95,-55},{95,-35}},
            lineColor={0,0,0},
            textString="level_start ="),
          Text(
            extent={{-95,50},{95,30}},
            lineColor={0,0,0},
            textString="level ="),
          Line(points={{-100,100},{-100,-100},{100,-100},{100,100}}, color={0,0,
                0})}),
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
<li><i>Dec. 12, 2008</i> by Ruediger Franke: replace energy and mass balances with
   common definition in BaseClasses.PartialLumpedVolume</li>
<li><i>Dec. 8, 2008</i> by Michael Wetter (LBNL):<br>
Implemented trace substances and missing equation for outflow of multi substance media at top port.</li>
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
      Diagram(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1},
          initialScale=0.2), graphics),
      uses(Modelica(version="2.2.1"), Modelica_Fluid(version="0.952")));
equation

    connect(heatPort, heatTransfer.heatPorts[1]) annotation (Line(
        points={{-100,0},{-87,0},{-87,8.88178e-016},{-74,8.88178e-016}},
        color={191,0,0},
        smooth=Smooth.None));
end Tank;

  package BaseClasses
    extends Modelica_Fluid.Icons.BaseClassLibrary;

      partial model PartialLumpedVolume
      "Lumped volume with dynamic mass and energy balance"
      import Modelica_Fluid.Types;
      import Modelica_Fluid.Types.Dynamics;

        outer Modelica_Fluid.System system "System properties";
        replaceable package Medium = 
          Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
            annotation (choicesAllMatching = true);

        // Assumptions
        parameter Types.Dynamics energyDynamics=system.energyDynamics
        "Formulation of energy balance" 
          annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
        parameter Types.Dynamics massDynamics=system.massDynamics
        "Formulation of mass balance" 
          annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
        final parameter Types.Dynamics substanceDynamics=massDynamics
        "Formulation of substance balance" 
          annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
        final parameter Types.Dynamics traceDynamics=massDynamics
        "Formulation of trace substance balance" 
          annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
        parameter Boolean use_d_nominal=energyDynamics<>Dynamics.SteadyState and massDynamics==Dynamics.SteadyState
        "= true if d_nominal is used for mass storage, else computed from medium"
          annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
        parameter Medium.Density d_nominal = Medium.density_pTX(Medium.p_default, Medium.T_default, Medium.X_default)
        "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
           annotation(Dialog(tab="Assumptions", group="Dynamics", enable=use_d_nominal));

        // Initialization
        parameter Medium.AbsolutePressure p_start = system.p_start
        "Start value of pressure" 
          annotation(Dialog(tab = "Initialization"));
        parameter Boolean use_T_start = true
        "= true, use T_start, otherwise h_start" 
          annotation(Dialog(tab = "Initialization"), Evaluate=true);
        parameter Medium.Temperature T_start=
          if use_T_start then system.T_start else Medium.temperature_phX(p_start,h_start,X_start)
        "Start value of temperature" 
          annotation(Dialog(tab = "Initialization", enable = use_T_start));
        parameter Medium.SpecificEnthalpy h_start=
          if use_T_start then Medium.specificEnthalpy_pTX(p_start, T_start, X_start) else Medium.h_default
        "Start value of specific enthalpy" 
          annotation(Dialog(tab = "Initialization", enable = not use_T_start));
        parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default
        "Start value of mass fractions m_i/m" 
          annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
        parameter Medium.ExtraProperty C_start[Medium.nC](
             quantity=Medium.extraPropertiesNames)=fill(0, Medium.nC)
        "Start value of trace substances" 
          annotation (Dialog(tab="Initialization", enable=Medium.nC > 0));

        Medium.BaseProperties medium(
          preferredMediumStates=true,
          p(start=p_start),
          h(start=h_start),
          T(start=T_start),
          Xi(start=X_start[1:Medium.nXi]));
        SI.Energy U "Internal energy of fluid";
        SI.Mass m "Mass of fluid";
        SI.Mass[Medium.nXi] mXi "Masses of independent components in the fluid";
        SI.Mass[Medium.nC] mC "Masses of trace substances in the fluid";
        // C need to be added here because unlike for Xi, which has medium.Xi,
        // there is no variable medium.C
        Medium.ExtraProperty C[Medium.nC] "Trace substance mixture content";

        // variables that need to be defined by an extending class
        SI.Volume fluidVolume "Volume";
        SI.MassFlowRate ms_flow "Mass flows across boundaries";
        SI.MassFlowRate[Medium.nXi] msXi_flow
        "Substance mass flows across boundaries";
        Medium.ExtraPropertyFlowRate[Medium.nC] msC_flow
        "Trace substance mass flows across boundaries";
        SI.EnthalpyFlowRate Hs_flow
        "Enthalpy flow across boundaries or energy source/sink";
        SI.HeatFlowRate Qs_flow
        "Heat flow across boundaries or energy source/sink";
        SI.Power Ws_flow "Work flow across boundaries or source term";
    protected
        parameter Boolean initialize_p = not Medium.singleState
        "= true to set up initial equations for pressure";
      equation

        // Total quantities
        if use_d_nominal then
          m = fluidVolume*d_nominal;
        else
          m = fluidVolume*medium.d;
        end if;
        mXi = m*medium.Xi;
        U = m*medium.u;
        mC = m*C;

        // Energy and mass balances
        if energyDynamics == Dynamics.SteadyState then
          0 = Hs_flow + Qs_flow + Ws_flow;
        else
          der(U) = Hs_flow + Qs_flow + Ws_flow;
        end if;

        if massDynamics == Dynamics.SteadyState then
          0 = ms_flow;
        else
          der(m) = ms_flow;
        end if;

        if substanceDynamics == Dynamics.SteadyState then
          zeros(Medium.nXi) = msXi_flow;
        else
          der(mXi) = msXi_flow;
        end if;

        if traceDynamics == Dynamics.SteadyState then
          zeros(Medium.nC)  = msC_flow;
        else
          der(mC)  = msC_flow;
        end if;

      initial equation
        // initialization of balances
        if energyDynamics == Dynamics.FixedInitial then
          if use_T_start then
            medium.T = T_start;
          else
            medium.h = h_start;
          end if;
        elseif energyDynamics == Dynamics.SteadyStateInitial then
          if use_T_start then
            der(medium.T) = 0;
          else
            der(medium.h) = 0;
          end if;
        end if;

        if massDynamics == Dynamics.FixedInitial then
          if initialize_p then
            medium.p = p_start;
          end if;
        elseif massDynamics == Dynamics.SteadyStateInitial then
          if initialize_p then
            der(medium.p) = 0;
          end if;
        end if;

        if substanceDynamics == Dynamics.FixedInitial then
          medium.Xi = X_start[1:Medium.nXi];
        elseif substanceDynamics == Dynamics.SteadyStateInitial then
          der(medium.Xi) = zeros(Medium.nXi);
        end if;

        if traceDynamics == Dynamics.FixedInitial then
          C = C_start[1:Medium.nC];
        elseif traceDynamics == Dynamics.SteadyStateInitial then
          der(C) = zeros(Medium.nC);
        end if;

        annotation (
          Documentation(info="<html>
Base class for an ideally mixed fluid volume with the ability to store mass and energy. 
The following source terms are part of the energy balance and must be specified in an extending class:
<ul>
<li><tt><b>Qs_flow</b></tt>, e.g. convective or latent heat flow rate across segment boundary, and</li> 
<li><tt><b>Ws_flow</b></tt>, work term, e.g. p*der(fluidVolume) if the volume is not constant.</li>
</ul>
The component volume <tt><b>fluidVolume</b></tt> is a variable which needs to be set in the extending class to complete the model.
<p>
Further source terms must be defined by an extending class for fluid flow across the segment boundary:
</p>
<ul>
<li><tt><b>Hs_flow</b></tt>, enthalpy flow,</li> 
<li><tt><b>ms_flow</b></tt>, mass flow,</li> 
<li><tt><b>msXi_flow</b></tt>, substance mass flow, and</li> 
<li><tt><b>msC_flow</b></tt>, trace substance mass flow.</li> 
</ul>
</html>"),Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
                -100},{100,100}}),
                  graphics));
      end PartialLumpedVolume;

      partial model PartialLumpedVolumePorts
      "Closed volume with a vector of fluid ports"
        extends PartialLumpedVolume;

      //Port definitions
        parameter Integer nPorts(min=1)=1 "Number of ports";
        Interfaces.FluidPorts_b[nPorts] ports(
                                      redeclare each package Medium = Medium)
        "Fluid outlets" 
          annotation (Placement(transformation(extent={{-10,-40},{10,40}},
            rotation=-90,
            origin={0,-100}),
            iconTransformation(extent={{-10,40},{10,-40}},
            rotation=-90,
            origin={0,-100})));
        Medium.AbsolutePressure ports_p_static
        "static pressure at the ports, inside the volume";

      //Transformation of kinetic energy
        parameter Boolean use_portDiameters=false
        "=true, kinetic energy and dissipation is accounted for in port pressure"
          annotation(Evaluate=true, Dialog(tab="Advanced",group="Ports"));
        parameter SI.Diameter portDiameters[nPorts] = fill(2.54e-2, nPorts)
        "Inner (hydraulic) diameters of ports (array)" 
          annotation(Dialog(tab="Advanced",group="Ports",enable= use_portDiameters));
        parameter Real[nPorts] zeta_in=fill(0, nPorts)
        "Hydraulic resistance into volume, 1 for total dissipation of kinetic energy and uniform flow distribution in pipe"
          annotation(Dialog(tab="Advanced",group="Ports",enable= use_portDiameters));
        parameter Real[nPorts] zeta_out=fill(1, nPorts)
        "Hydraulic resistance out of volume, 0 for ideal smooth outlet" 
          annotation(Dialog(tab="Advanced",group="Ports",enable= use_portDiameters));

        Medium.EnthalpyFlowRate ports_H_flow[nPorts];
        Medium.MassFlowRate ports_mXi_flow[nPorts,Medium.nXi];
        Medium.MassFlowRate[Medium.nXi] sum_ports_mXi_flow
        "Substance mass flows through ports";
        Medium.ExtraPropertyFlowRate ports_mC_flow[nPorts,Medium.nC];
        Medium.ExtraPropertyFlowRate[Medium.nC] sum_ports_mC_flow
        "Trace substance mass flows through ports";

        // Heat transfer through boundary
        parameter Boolean use_HeatTransfer = false
        "= true to use the HeatTransfer model" 
            annotation (Dialog(tab="Assumptions", group="Heat transfer"));
        replaceable model HeatTransfer = 
            Modelica_Fluid.Vessels.BaseClasses.HeatTransfer.IdealHeatTransfer 
          constrainedby
        Modelica_Fluid.Vessels.BaseClasses.HeatTransfer.PartialVesselHeatTransfer
        "Wall heat transfer" 
            annotation (Dialog(tab="Assumptions", group="Heat transfer",enable=use_HeatTransfer),choicesAllMatching=true);
        HeatTransfer heatTransfer(
          redeclare final package Medium = Medium,
          n=1,
          states = {medium.state}) 
            annotation (Placement(transformation(
              extent={{-10,-10},{30,30}},
              rotation=90,
              origin={-50,-10})));
        Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a heatPort if use_HeatTransfer 
          annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));

    protected
        parameter SI.Area[nPorts] portArea=Modelica.Constants.pi/4*{portDiameters[i]^2 for i in 1:nPorts};

      equation
        ms_flow = sum(ports.m_flow);
        msXi_flow = sum_ports_mXi_flow;
        msC_flow  = sum_ports_mC_flow;
        Hs_flow = sum(ports_H_flow);
        Qs_flow = heatTransfer.Q_flows[1];

        // Only one connection allowed to a port to avoid unwanted ideal mixing
        for i in 1:nPorts loop
          assert(cardinality(ports[i]) <= 1,"
each ports[i] of volume can at most be connected to one component.
If two or more connections are present, ideal mixing takes
place with these connections, which is usually not the intention
of the modeller. Increase nPorts to add an additional port.
");
        end for;

        // Boundary conditions
        for i in 1:nPorts loop
          if not use_portDiameters then
            ports[i].p = ports_p_static;
          else
            ports[i].p = ports_p_static - smooth(2, noEvent(ports[i].m_flow^2/(2*medium.d*
                  portArea[i]^2)*(if ports[i].m_flow < 0 then (1 + zeta_out[i]) else (1
                   - zeta_in[i]))));
          end if;
        end for;
        ports.h_outflow = fill(medium.h, nPorts);
        ports.Xi_outflow = fill(medium.Xi, nPorts);
        ports.C_outflow  = fill(C,         nPorts);

        for i in 1:nPorts loop
          ports_H_flow[i] = ports[i].m_flow * actualStream(ports[i].h_outflow)
          "Enthalpy flow";
          ports_mXi_flow[i,:] = ports[i].m_flow * actualStream(ports[i].Xi_outflow)
          "Component mass flow";
          ports_mC_flow[i,:]  = ports[i].m_flow * actualStream(ports[i].C_outflow)
          "Trace substance mass flow";
        end for;
        for i in 1:Medium.nXi loop
          sum_ports_mXi_flow[i] = sum(ports_mXi_flow[:,i]);
        end for;
        for i in 1:Medium.nC loop
          sum_ports_mC_flow[i]  = sum(ports_mC_flow[:,i]);
        end for;

        connect(heatPort, heatTransfer.heatPorts[1]) annotation (Line(
            points={{-100,0},{-87,0},{-87,8.88178e-016},{-74,8.88178e-016}},
            color={191,0,0},
            smooth=Smooth.None));
       annotation (
        Documentation(info="<html>
This base class extends PartialLumpedVolume by adding a vector of fluid ports 
and defining the respective source terms
<ul>
<li><tt>Hs_flow</tt>, enthalpy flow,</li> 
<li><tt>ms_flow</tt>, mass flow,</li> 
<li><tt>msXi_flow</tt>, substance mass flow, and</li> 
<li><tt>msC_flow</tt>, trace substance mass flow.</li> 
</ul>
An extending class still needs to define:
<ul>
<li><tt>Qs_flow</tt>, e.g. convective or latent heat flow rate across segment boundary,</li> 
<li><tt>Ws_flow</tt>, work term, e.g. p*der(V) if the volume is not constant, and</li>
<li><tt>V_lumped</tt>, the volume of the segment.</li>
</ul>
</html>"),
        Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},
                {100,100}}),
                graphics),
          Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},
                {100,100}}), graphics={Text(
              extent={{-150,110},{150,150}},
              textString="%name",
              lineColor={0,0,255})}));

      end PartialLumpedVolumePorts;

  partial model PartialDistributedVolume "Base class for a finite volume model"
      import Modelica_Fluid.Types;
      import Modelica_Fluid.Types.Dynamics;
    outer Modelica_Fluid.System system "System properties";

    replaceable package Medium = 
      Modelica.Media.Interfaces.PartialMedium "Medium in the component" 
        annotation (choicesAllMatching = true);

    // Discretization
    parameter Integer n=1 "Number of discrete flow volumes";

    // Assumptions
    parameter Types.Dynamics energyDynamics=system.energyDynamics
        "Formulation of energy balances" 
      annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
    parameter Types.Dynamics massDynamics=system.massDynamics
        "Formulation of mass balances" 
      annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
    final parameter Types.Dynamics substanceDynamics=massDynamics
        "Formulation of substance balances" 
      annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));
    final parameter Types.Dynamics traceDynamics=massDynamics
        "Formulation of trace substance balances" 
      annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));

    parameter Boolean use_d_nominal=energyDynamics<>Dynamics.SteadyState and massDynamics==Dynamics.SteadyState
        "= true if d_nominal is used for mass storage, else computed from medium"
      annotation(Evaluate=true, Dialog(tab = "Assumptions", group="Dynamics"));

    parameter Medium.Density d_nominal = Medium.density_pTX(Medium.p_default, Medium.T_default, Medium.X_default)
        "Nominal density (e.g. d_liquidWater = 995, d_air = 1.2)" 
       annotation(Dialog(tab="Assumptions", group="Dynamics", enable=use_d_nominal));

    //Initialization
    parameter Medium.AbsolutePressure p_a_start=system.p_start
        "Start value of pressure at port a" 
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.AbsolutePressure p_b_start=p_a_start
        "Start value of pressure at port b" 
      annotation(Dialog(tab = "Initialization"));
    final parameter Medium.AbsolutePressure[n] ps_start=if n > 1 then linspace(
          p_a_start, p_b_start, n) else {(p_a_start + p_b_start)/2}
        "Start value of pressure";

    parameter Boolean use_T_start=true "Use T_start if true, otherwise h_start"
       annotation(Evaluate=true, Dialog(tab = "Initialization"));
    parameter Medium.Temperature T_start=if use_T_start then system.T_start else 
                Medium.temperature_phX(
          (p_a_start + p_b_start)/2,
          h_start,
          X_start) "Start value of temperature" 
      annotation(Evaluate=true, Dialog(tab = "Initialization", enable = use_T_start));
    parameter Medium.SpecificEnthalpy h_start=if use_T_start then 
          Medium.specificEnthalpy_pTX(
          (p_a_start + p_b_start)/2,
          T_start,
          X_start) else Medium.h_default "Start value of specific enthalpy" 
      annotation(Evaluate=true, Dialog(tab = "Initialization", enable = not use_T_start));
    parameter Medium.MassFraction X_start[Medium.nX]=Medium.X_default
        "Start value of mass fractions m_i/m" 
      annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    parameter Medium.ExtraProperty C_start[Medium.nC](
         quantity=Medium.extraPropertiesNames)=fill(0, Medium.nC)
        "Start value of trace substances" 
      annotation (Dialog(tab="Initialization", enable=Medium.nC > 0));

    // Total quantities
    SI.Energy[n] Us "Internal energy of fluid";
    SI.Mass[n] ms "Fluid mass";
    SI.Mass[n,Medium.nXi] mXis "Substance mass";
    SI.Mass[n,Medium.nC] mCs "Trace substance mass";
    // C need to be added here because unlike for Xi, which has medium[:].Xi,
    // there is no variable medium[:].C
    Medium.ExtraProperty Cs[n, Medium.nC] "Trace substance mixture content";

    Medium.BaseProperties[n] mediums(
      each preferredMediumStates=true,
      p(start=ps_start),
      each h(start=h_start),
      each T(start=T_start),
      each Xi(start=X_start[1:Medium.nXi]));

    //Source terms, have to be set in inheriting class (to zero if not used)
    SI.Volume[n] fluidVolumes
        "Discretized volume, determine in inheriting class";
    Medium.MassFlowRate[n] ms_flows "Mass flow rate, source or sink";
    Medium.MassFlowRate[n,Medium.nXi] msXi_flows
        "Independent mass flow rates, source or sink";
    Medium.ExtraPropertyFlowRate[n,Medium.nC] msC_flows
        "Trace substance mass flow rates, source or sink";
    SI.EnthalpyFlowRate[n] Hs_flows "Enthalpy flow rate, source or sink";
    SI.HeatFlowRate[n] Qs_flows "Heat flow rate, source or sink";
    SI.Power[n] Ws_flows "Mechanical power, p*der(V) etc.";

    protected
    parameter Boolean initialize_p = not Medium.singleState
        "= true to set up initial equations for pressure";

  equation
    // Total quantities
    for i in 1:n loop
      if use_d_nominal then
        ms[i] =fluidVolumes[i]*d_nominal;
      else
        ms[i] =fluidVolumes[i]*mediums[i].d;
      end if;
      mXis[i, :] = ms[i]*mediums[i].Xi;
      mCs[i, :]  = ms[i]*Cs[i, :];
      Us[i] = ms[i]*mediums[i].u;
    end for;

    // Energy and mass balances
    if energyDynamics == Dynamics.SteadyState then
      for i in 1:n loop
        0 = Hs_flows[i] + Ws_flows[i] + Qs_flows[i];
      end for;
    else
      for i in 1:n loop
        der(Us[i]) = Hs_flows[i] + Ws_flows[i] + Qs_flows[i];
      end for;
    end if;
    if massDynamics == Dynamics.SteadyState then
      for i in 1:n loop
        0 = ms_flows[i];
      end for;
    else
      for i in 1:n loop
        der(ms[i]) = ms_flows[i];
      end for;
    end if;
    if substanceDynamics == Dynamics.SteadyState then
      for i in 1:n loop
        zeros(Medium.nXi) = msXi_flows[i, :];
      end for;
    else
      for i in 1:n loop
        der(mXis[i, :]) = msXi_flows[i, :];
      end for;
    end if;
    if traceDynamics == Dynamics.SteadyState then
      for i in 1:n loop
        zeros(Medium.nC)  = msC_flows[i, :];
      end for;
    else
      for i in 1:n loop
        der(mCs[i, :])  = msC_flows[i, :];
      end for;
    end if;

  initial equation
    // initialization of balances
    if energyDynamics == Dynamics.FixedInitial then
      if use_T_start then
        mediums.T = fill(T_start, n);
      else
        mediums.h = fill(h_start, n);
      end if;
    elseif energyDynamics == Dynamics.SteadyStateInitial then
      if use_T_start then
        der(mediums.T) = zeros(n);
      else
        der(mediums.h) = zeros(n);
      end if;
    end if;

    if massDynamics == Dynamics.FixedInitial then
      if initialize_p then
        mediums.p = ps_start;
      end if;
    elseif massDynamics == Dynamics.SteadyStateInitial then
      if initialize_p then
        der(mediums.p) = zeros(n);
      end if;
    end if;

    if substanceDynamics == Dynamics.FixedInitial then
      mediums.Xi = fill(X_start[1:Medium.nXi], n);
    elseif substanceDynamics == Dynamics.SteadyStateInitial then
      for i in 1:n loop
        der(mediums[i].Xi) = zeros(Medium.nXi);
      end for;
    end if;

    if traceDynamics == Dynamics.FixedInitial then
      Cs = fill(C_start[1:Medium.nC], n);
    elseif traceDynamics == Dynamics.SteadyStateInitial then
      for i in 1:n loop
        der(mCs[i,:])      = zeros(Medium.nC);
      end for;
    end if;

     annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
                -100},{100,100}}),
                         graphics),
                          Icon(coordinateSystem(preserveAspectRatio=true,
              extent={{-100,-100},{100,100}}), graphics),
        Documentation(info="<html>
Base class for <tt><b>n</b></tt> ideally mixed fluid volumes with the ability to store mass and energy.
It is inteded to model a one-dimensional spatial discretization of fluid flow according to the finite volume method. 
The following source terms are part of the energy balance and must be specified in an extending class:
<ul>
<li><tt><b>Qs_flows[n]</b></tt>, heat flow source term, e.g. conductive heat flows across segment boundaries, and</li> 
<li><tt><b>Ws_flows[n]</b></tt>, work source term.</li>
</ul>
The component volume <tt><b>fluidVolumes[n]</b></tt> is a variable which needs to be set in the extending class to complete the model.
<p>
Further source terms must be defined by an extending class for fluid flow across the segment boundary:
</p>
<ul>
<li><tt><b>Hs_flows[n]</b></tt>, enthalpy flow,</li> 
<li><tt><b>ms_flows[n]</b></tt>, mass flow,</li> 
<li><tt><b>msXi_flows[n]</b></tt>, substance mass flow, and</li> 
<li><tt><b>msC_flows[n]</b></tt>, trace substance mass flow.</li> 
</ul>
</html>"));
  end PartialDistributedVolume;

  package HeatTransfer

    partial model PartialVesselHeatTransfer
        "Base class for vessel heat transfer models"
      extends Modelica_Fluid.Interfaces.PartialHeatTransfer;

      input SI.Area[n] surfaceAreas "Heat transfer area";

      annotation(Documentation(info="<html>
Base class for vessel heat transfer models.
</html>"),Icon(coordinateSystem(preserveAspectRatio=true,  extent={{-100,-100},
                  {100,100}}), graphics={Ellipse(
                extent={{-60,64},{60,-56}},
                lineColor={0,0,0},
                fillPattern=FillPattern.Sphere,
                fillColor={232,0,0}), Text(
                extent={{-38,26},{40,-14}},
                lineColor={0,0,0},
                fillPattern=FillPattern.Sphere,
                fillColor={232,0,0},
                textString="%name")}));
    end PartialVesselHeatTransfer;

    model IdealHeatTransfer
        "IdealHeatTransfer: Ideal heat transfer without thermal resistance"
      extends PartialVesselHeatTransfer;
    equation
      Ts = heatPorts.T;
      annotation(Documentation(info="<html>
Ideal heat transfer without thermal resistance.
</html>"));
    end IdealHeatTransfer;

    model ConstantHeatTransfer
        "ConstantHeatTransfer: Constant heat transfer coefficient"
      extends PartialVesselHeatTransfer;
      parameter SI.CoefficientOfHeatTransfer alpha0
          "constant heat transfer coefficient";
      annotation(Documentation(info="<html>
Simple heat transfer correlation with constant heat transfer coefficient.
</html>"));
    equation
      Q_flows = {alpha0*surfaceAreas[i]*(heatPorts[i].T - Ts[i]) for i in 1:n};
    end ConstantHeatTransfer;
    annotation (Documentation(info="<html>
Heat transfer correlations for pipe models
</html>"));

  end HeatTransfer;

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
end Vessels;
