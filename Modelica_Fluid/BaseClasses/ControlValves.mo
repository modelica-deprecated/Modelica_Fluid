package ControlValves 
  partial model PartialValve "Base model for valves" 
    import Modelica_Fluid.Types.CvTypes;
    
  parameter Medium.AbsolutePressure pin_start = p_nom 
      "Start value of inlet pressure" 
    annotation(Dialog(tab = "Initialization"));
    
  extends BaseClasses.Common.PartialTwoPortTransport(
    medium_a(p(start=pin_start), T(start=T_start),
             h(start=h_start),   Xi(start=X_start[1:Medium.nXi])),
    medium_b(p(start=pout_start), T(start=T_start),
             h(start=h_start),   Xi(start=X_start[1:Medium.nXi])));
    
  parameter CvTypes.Temp CvData = CvTypes.Av "Selection of flow coefficient" 
     annotation(Dialog(group = "Flow Coefficient"));
  parameter SI.Area Av(fixed = if CvData==CvTypes.Av then true else false,
                       start = m_flow_nom/(sqrt(d_nom*dp_nom))*
                                           flowCharacteristic(stemPosition_nom)) = 0 
      "Av (metric) flow coefficient" 
     annotation(Dialog(group = "Flow Coefficient",
                       enable = (CvData==CvTypes.Av)));
  parameter Real Kv(unit="m3/h")=0 "Kv (metric) flow coefficient" 
    annotation(Dialog(group = "Flow Coefficient",
                      enable = (CvData==CvTypes.Kv)));
  parameter Real Cv(unit="USG/min")=0 "Cv (US) flow coefficient" 
    annotation(Dialog(group = "Flow Coefficient",
                      enable = (CvData==CvTypes.Cv)));
  parameter Medium.AbsolutePressure p_nom "Nominal inlet pressure" 
    annotation(Dialog(group="Nominal operating point"));
  parameter SI.Pressure dp_nom "Nominal pressure drop" 
    annotation(Dialog(group="Nominal operating point"));
  parameter Medium.MassFlowRate m_flow_nom "Nominal mass flowrate" 
    annotation(Dialog(group="Nominal operating point"));
  parameter Medium.Density d_nom = 1000 "Nominal inlet density" 
    annotation(Dialog(group="Nominal operating point"));
  parameter Real stemPosition_nom = 1 "Nominal stem position" 
    annotation(Dialog(group="Nominal operating point"));
  parameter Boolean CheckValve=false "Reverse flow stopped";
    
  replaceable function flowCharacteristic = 
      BaseClasses.ControlValves.ValveCharacteristics.linear 
    extends BaseClasses.ControlValves.ValveCharacteristics.baseFun 
      "Inherent flow characteristic" 
    annotation(choicesAllMatching=true);
  parameter Medium.AbsolutePressure pout_start = p_nom-dp_nom 
      "Start value of outlet pressure" 
    annotation(Dialog(tab = "Initialization"));
  parameter Boolean use_T_start = true "Use T_start if true, otherwise h_start"
    annotation(Dialog(tab = "Initialization"), Evaluate = true);
  parameter Medium.Temperature T_start=
    if use_T_start then Medium.T_default else Medium.temperature_phX(pin_start,h_start,X_start) 
      "Start value of inlet temperature" 
    annotation(Dialog(tab = "Initialization", enable = use_T_start));
  parameter Medium.SpecificEnthalpy h_start=
    if use_T_start then Medium.specificEnthalpy_pTX(pin_start, T_start, X_start[1:Medium.nXi]) else Medium.h_default 
      "Start value of specific enthalpy" 
    annotation(Dialog(tab = "Initialization", enable = not use_T_start));
  parameter Medium.MassFraction X_start[Medium.nX] = Medium.X_default 
      "Start value of mass fractions m_i/m" 
    annotation (Dialog(tab="Initialization", enable=Medium.nXi > 0));
    
  parameter Real delta=0.01 "Regularisation factor" annotation(Dialog(tab="Advanced"));
    
  Modelica.Blocks.Interfaces.RealInput stemPosition 
      "Stem position in the range 0-1" 
                                     annotation (extent=[-20,70; 20,110],   rotation=-90);
    
  Medium.Density d "Density at port a";
  Medium.Temperature T "Temperature at port a";
  protected 
  function sqrtR = Utilities.regRoot(delta = delta*dp_nom);
  annotation (
    Icon(Text(extent=[-143,-66; 148,-106],  string="%name"),
      Line(points=[0,60; 0,0],   style(
          color=0,
          thickness=2,
          fillPattern=1)),
      Polygon(points=[-100,50; -100,-50; 0,0; -100,50],  style(
          color=0,
          thickness=2,
          fillPattern=1)),
      Polygon(points=[100,50; 0,0; 100,-50; 100,50],  style(
          color=0,
          thickness=2,
          fillPattern=1)),
      Rectangle(extent=[-20,70; 20,50],   style(
          color=0,
          fillColor=0,
          fillPattern=1))),
    Diagram,
    Documentation(info="<HTML>
<p>This is the base model for the <tt>ValveIncompressible</tt>, <tt>ValveVaporizing</tt>, and <tt>ValveCompressible</tt> valve models. The model is based on the IEC 534 / ISA S.75 standards for valve sizing.
<p>The model optionally supports reverse flow conditions (assuming symmetrical behaviour) or check valve operation, and has been suitably modified to avoid numerical singularities at zero pressure drop. 
<p><b>Modelling options</b></p>
<p>The following options are available to specify the valve flow coefficient in fully open conditions:
<ul><li><tt>CvData = Modelica_Fluid.Types.CvTypes.Av</tt>: the flow coefficient is given by the metric <tt>Av</tt> coefficient (m^2).
<li><tt>CvData = Modelica_Fluid.Types.CvTypes.Kv</tt>: the flow coefficient is given by the metric <tt>Kv</tt> coefficient (m^3/h).
<li><tt>CvData = Modelica_Fluid.Types.CvTypes.Cv</tt>: the flow coefficient is given by the US <tt>Cv</tt> coefficient (USG/min).
<li><tt>CvData = Modelica_Fluid.Types.CvTypes.OpPoint</tt>: the flow is computed from the nominal operating point specified by <tt>p_nom</tt>, <tt>dp_nom</tt>, <tt>m_flow_nom</tt>, <tt>d_nom</tt>, <tt>stemPosition_nom</tt>.
</ul>
<p>The nominal pressure drop <tt>dp_nom</tt> must always be specified; to avoid numerical singularities, the flow characteristic is modified for pressure drops less than <tt>b*dp_nom</tt> (the default value is 1% of the nominal pressure drop). Increase this parameter if numerical problems occur in valves with very low pressure drops.
<p>If <tt>CheckValve</tt> is true, then the flow is stopped when the outlet pressure is higher than the inlet pressure; otherwise, reverse flow takes place.
<p>The inherent flow characteristic <tt>flowCharacteristic</tt>, linear by default, can be replaced by any user-defined function (e.g. equal percentage, quick opening, etc.).
</HTML>",
      revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"),
    Coordsys(grid=[2,2], scale=0));
  initial equation 
  if CvData == CvTypes.Kv then
    Av = 2.7778e-5*Kv "Unit conversion";
  elseif CvData == CvTypes.Cv then
    Av = 2.4027e-5*Cv "Unit conversion";
  end if;
  assert(CvData>=0 and CvData<=3, "Invalid CvData");
  equation 
  T = medium_a.T;
  d = medium_a.d;
  end PartialValve;
  
package ValveCharacteristics "Functions for valve characteristics" 
  partial function baseFun "Base class for valve characteristics" 
    extends Modelica.Icons.Function;
    input Real pos "Stem position (per unit)";
    output Real rc "Relative flow coefficient (per unit)";
  end baseFun;
    
  function linear "Linear characteristic" 
    extends baseFun;
  algorithm 
    rc := pos;
  end linear;
    
  function one "Constant characteristic" 
    extends baseFun;
  algorithm 
    rc := 1;
  end one;
    
  function quadratic "Quadratic characteristic" 
    extends baseFun;
  algorithm 
    rc := pos*pos;
  end quadratic;
    
  function equalPercentage "Equal percentage characteristic" 
    extends baseFun;
    input Real rangeability = 20 "Rangeability";
    input Real delta = 0.01;
  algorithm 
    rc := if pos > delta then rangeability^(pos-1) else 
            pos/delta*rangeability^(delta-1);
    annotation (Documentation(info="<html>
This characteristic is such that the relative change of the flow coefficient is proportional to the change in the stem position:
<p> d(rc)/d(pos) = k d(pos).
<p> The constant k is expressed in terms of the rangeability, i.e. the ratio between the maximum and the minimum useful flow coefficient:
<p> rangeability = exp(k) = rc(1.0)/rc(0.0).
<p> The theoretical characteristic has a non-zero opening when pos = 0; the implemented characteristic is modified so that the valve closes linearly when pos &lt delta.
</html>"));
  end equalPercentage;
    
end ValveCharacteristics;
end ControlValves;
