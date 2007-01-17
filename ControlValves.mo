package ControlValves "Various variants of valve components" 
    extends Modelica_Fluid.Icons.VariantLibrary;
  
    model ValveIncompressible "Valve for (almost) incompressible fluids" 
    extends BaseClasses.PartialValve;
    import Modelica_Fluid.Types.CvTypes;
    annotation (
    Icon,
    Diagram,
    Documentation(info="<HTML>
<p>Valve model according to the IEC 534/ISA S.75 standards for valve sizing, incompressible fluids. <p>
Extends the <tt>BaseClasses.ControlValves.PartialValve</tt> model (see the corresponding documentation for common valve features).
<p>This model can be used with any low compressibility fluids, such as liquids or gases at very low pressure drops.</p>
<p>If <tt>CheckValve</tt> is false, the valve supports reverse flow, with a symmetric flow characteric curve. Otherwise, reverse flow is stopped (check valve behaviour).</p>
</html>",
      revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"),
      Coordsys(grid=[1,1], scale=0));
    initial equation 
    if CvData == CvTypes.OpPoint then
      m_flow_nom = flowCharacteristic(stemPosition_nom)*Av*sqrt(d_nom)*sqrtR(dp_nom) 
        "Determination of Av by the operating point";
    end if;
    
    equation 
    if CheckValve then
      m_flow = flowCharacteristic(stemPosition)*Av*sqrt(d)*
               smooth(0,if dp>=0 then sqrtR(dp) else 0);
    else
      m_flow = flowCharacteristic(stemPosition)*Av*sqrt(d)*sqrtR(dp);
    end if;
    end ValveIncompressible;
  
  model ValveVaporizing 
    "Valve for possibly vaporizing (almost) incompressible fluids, accounts for choked flow conditions" 
    import Modelica_Fluid.Types.CvTypes;
    extends BaseClasses.PartialValve(
      redeclare replaceable package Medium = 
      Modelica.Media.Interfaces.PartialTwoPhaseMedium);
    parameter Real Fl_nom=0.9 "Liquid pressure recovery factor";
    replaceable function FlCharacteristic = 
        Modelica_Fluid.ControlValves.BaseClasses.ValveCharacteristics.one 
      extends 
      Modelica_Fluid.ControlValves.BaseClasses.ValveCharacteristics.baseFun 
      "Pressure recovery characteristic";
    Real Ff "Ff coefficient (see IEC/ISA standard)";
    Real Fl "Pressure recovery coefficient Fl (see IEC/ISA standard)";
    Medium.AbsolutePressure pv "Saturation pressure";
    SI.Pressure dpEff "Effective pressure drop";
    Medium.AbsolutePressure pin "Inlet pressure";
    Medium.AbsolutePressure pout "Outlet pressure";
    annotation (
      Icon,
      Diagram,
      Documentation(info="<HTML>
<p>Valve model according to the IEC 534/ISA S.75 standards for valve sizing, incompressible fluid at the inlet, and possibly two-phase fluid at the outlet, with resulting choked flow conditions. <p>
Extends the <tt>BaseClasses.ControlValves.PartialValve</tt> model (see the corresponding documentation for common valve features).<p>
The model operating range includes choked flow operation, which takes place for low outlet pressures due to flashing in the vena contracta; otherwise, non-choking conditions are assumed.
<p>This model must be used with two-phase medium models, to describe the liquid and (possible) two-phase conditions.
<p>The default liquid pressure recovery coefficient <tt>Fl</tt> is constant and given by the parameter <tt>Fl_nom</tt>. The relative change (per unit) of the recovery coefficient can be specified as a given function of the valve opening by replacing the <tt>FlCharacteristic</tt> function.
<p>If <tt>CheckValve</tt> is false, the valve supports reverse flow, with a symmetric flow characteric curve. Otherwise, reverse flow is stopped (check valve behaviour).</p>

</HTML>",
        revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"));
  initial equation 
    if CvData == CvTypes.OpPoint then
      m_flow_nom = flowCharacteristic(stemPosition_nom)*Av*sqrt(d_nom)*sqrtR(dp_nom) 
        "Determination of Av by the operating point";
    end if;
  equation 
    pin = port_a.p;
    pout = port_b.p;
    pv = Medium.saturationPressure_sat(medium_a.sat);
    Ff = 0.96 - 0.28*sqrt(pv/Medium.fluidConstants[1].criticalPressure);
    Fl = Fl_nom*FlCharacteristic(stemPosition);
    dpEff = if pout < (1 - Fl^2)*pin + Ff*Fl^2*pv then 
              Fl^2*(pin - Ff*pv) else dp 
      "Effective pressure drop, accounting for possible choked conditions";
    if CheckValve then
       m_flow = flowCharacteristic(stemPosition)*Av*sqrt(d)*
           (if dpEff>=0 then sqrtR(dpEff) else 0);
     else
       m_flow = flowCharacteristic(stemPosition)*Av*sqrt(d)*sqrtR(dpEff);
    end if;
  end ValveVaporizing;
  
  model ValveCompressible 
    "Valve for compressible fluids, accounts for choked flow conditions" 
    extends BaseClasses.PartialValve;
    import Modelica_Fluid.Types.CvTypes;
    parameter Real Fxt_full=0.5 "Fk*xt critical ratio at full opening";
    replaceable function xtCharacteristic = 
        Modelica_Fluid.ControlValves.BaseClasses.ValveCharacteristics.one 
      extends 
      Modelica_Fluid.ControlValves.BaseClasses.ValveCharacteristics.baseFun 
      "Critical ratio characteristic";
    Real Fxt;
    Real x "Pressure drop ratio";
    Real xs "Saturated pressure drop ratio";
    Real Y "Compressibility factor";
    Medium.AbsolutePressure p "Inlet pressure";
  protected 
    parameter Real Fxt_nom(fixed=false) "Nominal Fxt";
    parameter Real x_nom(fixed=false) "Nominal pressure drop ratio";
    parameter Real xs_nom(fixed=false) "Nominal saturated pressure drop ratio";
    parameter Real Y_nom(fixed=false) "Nominal compressibility factor";
    
    annotation (
    Icon,
    Diagram,
    Documentation(info="<HTML>
<p>Valve model according to the IEC 534/ISA S.75 standards for valve sizing, compressible fluid, no phase change, including choked conditions. <p>
Extends the <tt>BaseClasses.ControlValves.PartialValve</tt> model (see the corresponding documentation for common valve features).
<p>This model can be used with gases at moderate to high pressure ratios.</p>

<p>The product Fk*xt is given by the parameter <tt>Fxt_full</tt>, and is assumed constant by default. The relative change (per unit) of the xt coefficient with the valve opening can be specified by replacing the <tt>xtCharacteristic</tt> function.
<p>If <tt>CheckValve</tt> is false, the valve supports reverse flow, with a symmetric flow characteric curve. Otherwise, reverse flow is stopped (check valve behaviour).</p>

</HTML>",
      revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"));
  initial equation 
    if CvData == CvTypes.OpPoint then
      // Determination of Av by the nominal operating point conditions
      Fxt_nom = Fxt_full*xtCharacteristic(stemPosition_nom);
      x_nom = dp_nom/p_nom;
      xs_nom = smooth(0, if x_nom > Fxt_nom then Fxt_nom else x_nom);
      Y_nom = 1 - abs(xs_nom)/(3*Fxt_nom);
      m_flow_nom = flowCharacteristic(stemPosition_nom)*Av*Y_nom*sqrt(d_nom)*sqrtR(p_nom*xs_nom);
    else
      // Dummy values
      Fxt_nom = 0;
      x_nom = 0;
      xs_nom = 0;
      Y_nom = 0;
    end if;
    
  equation 
    p = noEvent(if dp>=0 then port_a.p else port_b.p);
    Fxt = Fxt_full*xtCharacteristic(stemPosition);
    x = dp/p;
    xs = smooth(0, if x < -Fxt then -Fxt else if x > Fxt then Fxt else x);
    Y = 1 - abs(xs)/(3*Fxt);
    if CheckValve then
      m_flow = flowCharacteristic(stemPosition)*Av*Y*sqrt(d)*
        smooth(0,if xs>=0 then sqrtR(p*xs) else 0);
    else
      m_flow = flowCharacteristic(stemPosition)*Av*Y*sqrt(d)*sqrtR(p*xs);
    end if;
  end ValveCompressible;
  
  model ValveLinear "Valve for water/steam flows with linear pressure drop" 
    extends Modelica_Fluid.Interfaces.PartialTwoPortTransport;
    parameter Types.HydraulicConductance Kv 
      "Hydraulic conductance at full opening";
    Modelica.Blocks.Interfaces.RealInput opening 
    annotation (extent=[-20,70; 20,110],   rotation=-90);
  equation 
    m_flow = Kv*opening*dp;
    
  annotation (
    Icon(
        Polygon(points=[-100,50; -100,-50; 0,0; -100,50],  style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Line(points=[0,60; 0,0],   style(
            color=0,
            thickness=2,
            fillPattern=1)),
        Rectangle(extent=[-20,70; 20,50],   style(
            color=0,
            fillColor=0,
            fillPattern=1)),
        Polygon(points=[100,50; 0,0; 100,-50; 100,50],  style(
            color=0,
            thickness=2,
            fillPattern=1)),
           Text(extent=[-143,-66; 148,-106],  string="%name")),
    Diagram,
    Documentation(info="<HTML>
<p>This very simple model provides a pressure drop which is proportional to the flowrate and to the <tt>opening</tt> signal, without computing any fluid property.
<p>A medium model must be nevertheless be specified, so that the fluid ports can be connected to other components using the same medium model.
</HTML>",
      revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"),
      Coordsys(grid=[1,1], scale=0));
  end ValveLinear;
  
  model ValveDiscrete "Valve for water/steam flows with linear pressure drop" 
    extends Modelica_Fluid.Interfaces.PartialTwoPortTransport;
    parameter Modelica_Fluid.Types.HydraulicConductance Kv 
      "Hydraulic conductance for open valve (m_flow = Kv*dp)";
    parameter Real Kv_small_rel = 0 
      "Relative hydraulic conductance for closed valve (m_flow = Kv_small_rel*Kv*dp)";
    Modelica.Blocks.Interfaces.BooleanInput open 
    annotation (extent=[-20,60; 20,100],   rotation=-90);
  equation 
    m_flow = if open then Kv*dp else Kv_small_rel*Kv*dp;
    
  annotation (
    Icon(
        Line(points=[0,50; 0,0],   style(
            color=0,
            rgbcolor={0,0,0},
            fillPattern=1)),
        Rectangle(extent=[-20,60; 20,50],   style(
            color=0,
            fillColor=0,
            fillPattern=1)),
           Text(extent=[-145,-58; 146,-98],   string="%name"),
        Polygon(points=[-100,50; 100,-50; 100,50; 0,0; -100,-50; -100,50], style(
            color=0,
            rgbcolor={0,0,0},
            fillColor=DynamicSelect(7, if open > 0.5 then 2 else 7)))),
    Diagram,
    Documentation(info="<HTML>
<
<p>
This very simple model provides a pressure drop which is proportional to the flowrate if the Boolean open signal is <b>true</b>. Otherwise, the
mass flow rate is zero. If Kv_small_rel > 0, a small leakage
mass flow rate occurs when open = <b>false</b>. This might be
useful in certain situations when the model is not
mathematically well-defined due to a closed valve.
</p>
<p>
In a diagram animation, the valve is shown in \"green\", when
it is open.
</p>
</HTML>",
      revisions="<html>
<ul>
<li><i>Nov 2005</i>
    by Katja Poschlad (based on ValveLinear).</li>
</ul>
</html>"),
      Coordsys(grid=[1,1], scale=0));
  end ValveDiscrete;
  
  package BaseClasses 
    extends Modelica_Fluid.Icons.BaseClassLibrary;
    partial model PartialValve "Base model for valves" 
      import Modelica_Fluid.Types.CvTypes;
      import Modelica_Fluid;
      
    parameter Medium.AbsolutePressure pin_start = p_nom 
        "Start value of inlet pressure" 
      annotation(Dialog(tab = "Initialization"));
      
    extends Modelica_Fluid.Interfaces.PartialTwoPortTransport(medium_a(
          p(start=pin_start),
          T(start=T_start),
          h(start=h_start),
          Xi(start=X_start[1:Medium.nXi])), medium_b(
          p(start=pout_start),
          T(start=T_start),
          h(start=h_start),
          Xi(start=X_start[1:Medium.nXi])),
          m_flow(start = m_flow_start));
      
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
        Modelica_Fluid.ControlValves.BaseClasses.ValveCharacteristics.linear 
        extends 
        Modelica_Fluid.ControlValves.BaseClasses.ValveCharacteristics.baseFun 
        "Inherent flow characteristic" 
      annotation(choicesAllMatching=true);
    parameter Medium.AbsolutePressure pout_start = p_nom-dp_nom 
        "Start value of outlet pressure" 
      annotation(Dialog(tab = "Initialization"));
    parameter Medium.MassFlowRate m_flow_start = m_flow_nom 
        "Start value of mass flow rate" 
      annotation(Dialog(tab = "Initialization"));
    parameter Boolean use_T_start = true 
        "Use T_start if true, otherwise h_start" 
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
  end BaseClasses;
  annotation (Documentation(info="<html>
 
</html>"));
end ControlValves;
