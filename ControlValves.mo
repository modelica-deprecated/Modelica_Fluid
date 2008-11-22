within Modelica_Fluid;
package ControlValves "Various variants of valve components"
    extends Modelica_Fluid.Icons.VariantLibrary;

    model ValveIncompressible "Valve for (almost) incompressible fluids"
      extends BaseClasses.PartialValve;
    import Modelica_Fluid.Types.CvTypes;
    annotation (
    Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics),
    Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics),
    Documentation(info="<HTML>
<p>Valve model according to the IEC 534/ISA S.75 standards for valve sizing, incompressible fluids. <p>
Extends the <tt>BaseClasses.ControlValves.PartialValve</tt> model (see the corresponding documentation for common valve features).
<p>This model can be used with any low compressibility fluids, such as liquids or gases at very low pressure drops.</p>
<p>If <tt>CheckValve</tt> is false, the valve supports reverse flow, with a symmetric flow characteric curve. Otherwise, reverse flow is stopped (check valve behaviour).</p>
 
<p>
The treatment of parameters <b>Kv</b> and <b>Cv</b> is
explained in detail in the 
<a href=\"Modelica://Modelica_Fluid.UsersGuide.ComponentDefinition.ValveCharacteristics\">Users Guide</a>.
</p>
 
</html>",
      revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"));
    initial equation
      if CvData == CvTypes.OpPoint then
          m_flow_nominal = valveCharacteristic(stemPosition_nominal)*Av*sqrt(d_nominal)*Utilities.regRoot(dp_nominal, delta*dp_nominal)
        "Determination of Av by the operating point";
      end if;

    equation
      if CheckValve then
          m_flow = valveCharacteristic(modifiedStemPosition)*Av*sqrt(port_a_d_inflow)*
                      smooth(0,if dp>=0 then Utilities.regRoot(dp, delta*dp_nominal) else 0);
      else
        // m_flow = valveCharacteristic(stemPosition)*Av*sqrt(d)*sqrtR(dp);
        m_flow = valveCharacteristic(modifiedStemPosition)*Av*sqrt(port_a_d_inflow)*
          Utilities.regRoot(dp, delta*dp_nominal);
      end if;
    end ValveIncompressible;

  model ValveVaporizing
    "Valve for possibly vaporizing (almost) incompressible fluids, accounts for choked flow conditions"
    import Modelica_Fluid.Types.CvTypes;
    extends BaseClasses.PartialValve(
      final compute_T = true,
      redeclare replaceable package Medium = 
          Modelica.Media.Water.WaterIF97_ph                                    constrainedby
        Modelica.Media.Interfaces.PartialTwoPhaseMedium);
    parameter Real Fl_nominal=0.9 "Liquid pressure recovery factor";
    replaceable function FlCharacteristic = 
        Modelica_Fluid.ControlValves.BaseClasses.ValveCharacteristics.one 
      constrainedby
      Modelica_Fluid.ControlValves.BaseClasses.ValveCharacteristics.baseFun
      "Pressure recovery characteristic";
    Real Ff "Ff coefficient (see IEC/ISA standard)";
    Real Fl "Pressure recovery coefficient Fl (see IEC/ISA standard)";
    Medium.AbsolutePressure pv "Saturation pressure";
    SI.Pressure dpEff "Effective pressure drop";
    Medium.AbsolutePressure pin "Inlet pressure";
    Medium.AbsolutePressure pout "Outlet pressure";
    annotation (
      Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}),
           graphics),
      Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}),
              graphics),
      Documentation(info="<HTML>
<p>Valve model according to the IEC 534/ISA S.75 standards for valve sizing, incompressible fluid at the inlet, and possibly two-phase fluid at the outlet, with resulting choked flow conditions. <p>
Extends the <tt>BaseClasses.ControlValves.PartialValve</tt> model (see the corresponding documentation for common valve features).<p>
The model operating range includes choked flow operation, which takes place for low outlet pressures due to flashing in the vena contracta; otherwise, non-choking conditions are assumed.
<p>This model must be used with two-phase medium models, to describe the liquid and (possible) two-phase conditions.
<p>The default liquid pressure recovery coefficient <tt>Fl</tt> is constant and given by the parameter <tt>Fl_nominal</tt>. The relative change (per unit) of the recovery coefficient can be specified as a given function of the valve opening by replacing the <tt>FlCharacteristic</tt> function.
<p>If <tt>CheckValve</tt> is false, the valve supports reverse flow, with a symmetric flow characteric curve. Otherwise, reverse flow is stopped (check valve behaviour).</p>
 
<p>
The treatment of parameters <b>Kv</b> and <b>Cv</b> is
explained in detail in the 
<a href=\"Modelica://Modelica_Fluid.UsersGuide.ComponentDefinition.ValveCharacteristics\">Users Guide</a>.
</p>
 
</HTML>",
        revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"));
  initial equation
    assert(not CvData == CvTypes.OpPoint, "OpPoint option not supported for vaporizing valve");
  equation
    pin = port_a.p;
    pout = port_b.p;
    pv = Medium.saturationPressure(port_a_T);
    Ff = 0.96 - 0.28*sqrt(pv/Medium.fluidConstants[1].criticalPressure);
    Fl = Fl_nominal*FlCharacteristic(stemPosition);
    dpEff = if pout < (1 - Fl^2)*pin + Ff*Fl^2*pv then 
              Fl^2*(pin - Ff*pv) else dp
      "Effective pressure drop, accounting for possible choked conditions";
    if CheckValve then
       m_flow = valveCharacteristic(modifiedStemPosition)*Av*sqrt(port_a_d_inflow)*
           smooth(0,if dpEff>=0 then sqrtR(dpEff) else 0);
     else
       // m_flow = valveCharacteristic(stemPosition)*Av*sqrt(d)*sqrtR(dpEff);
       m_flow = valveCharacteristic(modifiedStemPosition)*Av*sqrt(port_a_d_inflow)*sqrtR(dpEff);
    end if;
    assert(m_flow > -0.1 * m_flow_nominal, "Too big backflow");
  end ValveVaporizing;

  model ValveCompressible
    "Valve for compressible fluids, accounts for choked flow conditions"
    extends BaseClasses.PartialValve;
    import Modelica_Fluid.Types.CvTypes;
    parameter Medium.AbsolutePressure p_nominal "Nominal inlet pressure" 
    annotation(Dialog(group="Nominal operating point"));
    parameter Real Fxt_full=0.5 "Fk*xt critical ratio at full opening";
    replaceable function xtCharacteristic = 
        Modelica_Fluid.ControlValves.BaseClasses.ValveCharacteristics.one 
      constrainedby
      Modelica_Fluid.ControlValves.BaseClasses.ValveCharacteristics.baseFun
      "Critical ratio characteristic";
    Real Fxt;
    Real x "Pressure drop ratio";
    Real xs "Saturated pressure drop ratio";
    Real Y "Compressibility factor";
    Medium.AbsolutePressure p "Inlet pressure";
  protected
    parameter Real Fxt_nominal(fixed=false) "Nominal Fxt";
    parameter Real x_nominal(fixed=false) "Nominal pressure drop ratio";
    parameter Real xs_nominal(fixed=false)
      "Nominal saturated pressure drop ratio";
    parameter Real Y_nominal(fixed=false) "Nominal compressibility factor";

    annotation (
    Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,
              100}}),
         graphics),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
              100,100}}),
            graphics),
    Documentation(info="<HTML>
<p>Valve model according to the IEC 534/ISA S.75 standards for valve sizing, compressible fluid, no phase change, including choked conditions. <p>
Extends the <tt>BaseClasses.ControlValves.PartialValve</tt> model (see the corresponding documentation for common valve features).
<p>This model can be used with gases at moderate to high pressure ratios.</p>
 
<p>The product Fk*xt is given by the parameter <tt>Fxt_full</tt>, and is assumed constant by default. The relative change (per unit) of the xt coefficient with the valve opening can be specified by replacing the <tt>xtCharacteristic</tt> function.
<p>If <tt>CheckValve</tt> is false, the valve supports reverse flow, with a symmetric flow characteric curve. Otherwise, reverse flow is stopped (check valve behaviour).</p>
 
<p>
The treatment of parameters <b>Kv</b> and <b>Cv</b> is
explained in detail in the 
<a href=\"Modelica://Modelica_Fluid.UsersGuide.ComponentDefinition.ValveCharacteristics\">Users Guide</a>.
</p>
 
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
      Fxt_nominal = Fxt_full*xtCharacteristic(stemPosition_nominal);
      x_nominal = dp_nominal/p_nominal;
      xs_nominal = smooth(0, if x_nominal > Fxt_nominal then Fxt_nominal else x_nominal);
      Y_nominal = 1 - abs(xs_nominal)/(3*Fxt_nominal);
      m_flow_nominal = valveCharacteristic(stemPosition_nominal)*Av*Y_nominal*sqrt(d_nominal)*sqrtR(p_nominal*xs_nominal);
    else
      // Dummy values
      Fxt_nominal = 0;
      x_nominal = 0;
      xs_nominal = 0;
      Y_nominal = 0;
    end if;

  equation
    p = noEvent(if dp>=0 then port_a.p else port_b.p);
    Fxt = Fxt_full*xtCharacteristic(modifiedStemPosition);
    x = dp/p;
    xs = smooth(0, if x < -Fxt then -Fxt else if x > Fxt then Fxt else x);
    Y = 1 - abs(xs)/(3*Fxt);
    if CheckValve then
      m_flow = valveCharacteristic(modifiedStemPosition)*Av*Y*sqrt(port_a_d_inflow)*
        smooth(0,if xs>=0 then sqrtR(p*xs) else 0);
    else
      // m_flow = valveCharacteristic(stemPosition)*Av*Y*sqrt(d)*sqrtR(p*xs);
      m_flow = valveCharacteristic(modifiedStemPosition)*Av*Y*
                    Modelica_Fluid.Utilities.regRoot2(p*xs, delta*dp_nominal, port_a_d_inflow, port_b_d_inflow);
    end if;
  end ValveCompressible;

  model ValveLinear "Valve for water/steam flows with linear pressure drop"
    extends Modelica_Fluid.PressureLosses.BaseClasses.PartialTwoPortTransport;
    parameter SI.AbsolutePressure dp_nominal
      "Nominal pressure drop at full opening";
    parameter Medium.MassFlowRate m_flow_nominal
      "Nominal mass flowrate at full opening";
    final parameter Types.HydraulicConductance Kv = m_flow_nominal/dp_nominal
      "Hydraulic conductance at full opening";
    Modelica.Blocks.Interfaces.RealInput opening(min=0,max=1)
      "=1: completely open, =0: completely closed" 
    annotation (Placement(transformation(
          origin={0,90},
          extent={{-20,-20},{20,20}},
          rotation=270), iconTransformation(
          extent={{-20,-20},{20,20}},
          rotation=270,
          origin={0,80})));
    parameter Real minOpening(min=0, max=0.1)=0
      "Minimum position of opening (leakage flow to improve numerics)" 
    annotation(Dialog(tab="Advanced"));
    Real modifiedOpening
      "Modified, actually used opening, so that the valve is not completely closed to improve numerics";

  equation
    modifiedOpening = smooth(0,noEvent(if opening > minOpening then opening else minOpening));
    m_flow = Kv*modifiedOpening*dp;

  annotation (
    Icon(coordinateSystem(
          preserveAspectRatio=true,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics={
          Polygon(
            points={{-100,50},{-100,-50},{0,0},{-100,50}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Line(points={{0,50},{0,0}}, color={0,0,0}),
          Rectangle(
            extent={{-20,60},{20,50}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{100,50},{0,0},{100,-50},{100,50}},
            lineColor={0,0,0},
            fillColor={255,255,255},
            fillPattern=FillPattern.Solid),
          Line(
            points={{30,-50},{-60,-50}},
            color={0,128,255},
            smooth=Smooth.None),
          Polygon(
            points={{20,-36},{60,-50},{20,-64},{20,-36}},
            lineColor={0,128,255},
            smooth=Smooth.None,
            fillColor={0,128,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-150,-80},{150,-120}},
            textString="%name",
            lineColor={0,0,255})}),
    Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={1,1}), graphics),
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
</html>"));
  end ValveLinear;

  model ValveDiscrete "Valve for water/steam flows with linear pressure drop"
    extends Modelica_Fluid.PressureLosses.BaseClasses.PartialTwoPortTransport;
    parameter SI.Pressure dp_nominal "Nominal pressure drop at full opening";
    parameter Medium.MassFlowRate m_flow_nominal
      "Nominal mass flowrate at full opening";
    final parameter Types.HydraulicConductance Kv = m_flow_nominal/dp_nominal
      "Hydraulic conductance at full opening";
    parameter Real Kv_small_rel(min=0, max=0.1) = 0
      "Relative hydraulic conductance for closed valve (m_flow = Kv_small_rel*Kv*dp)";
    Modelica.Blocks.Interfaces.BooleanInput open 
    annotation (Placement(transformation(
          origin={0,80},
          extent={{-20,-20},{20,20}},
          rotation=270)));
  equation
    m_flow = if open then Kv*dp else Kv_small_rel*Kv*dp;

  annotation (
    Icon(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={2,2}), graphics={
          Line(points={{0,50},{0,0}}, color={0,0,0}),
          Rectangle(
            extent={{-20,60},{20,50}},
            lineColor={0,0,0},
            fillColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Polygon(
            points={{-100,50},{100,-50},{100,50},{0,0},{-100,-50},{-100,50}},
            fillColor={255,255,255},
            lineColor={0,0,0},
            fillPattern=FillPattern.Solid),
          Line(
            points={{30,-50},{-60,-50}},
            color={0,128,255},
            smooth=Smooth.None),
          Polygon(
            points={{20,-36},{60,-50},{20,-64},{20,-36}},
            lineColor={0,128,255},
            smooth=Smooth.None,
            fillColor={0,128,255},
            fillPattern=FillPattern.Solid),
          Text(
            extent={{-150,-80},{150,-120}},
            textString="%name",
            lineColor={0,0,255})}),
    Diagram(coordinateSystem(
          preserveAspectRatio=false,
          extent={{-100,-100},{100,100}},
          grid={2,2}), graphics),
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
</html>"));

  end ValveDiscrete;

  package BaseClasses
    extends Modelica_Fluid.Icons.BaseClassLibrary;
    partial model PartialValve "Base model for valves"

      import Modelica_Fluid.Types.CvTypes;
      extends Modelica_Fluid.PressureLosses.BaseClasses.PartialTwoPortTransport(
          dp_start=dp_nominal, m_flow_start=m_flow_nominal);
      parameter CvTypes CvData=CvTypes.OpPoint "Selection of flow coefficient" 
       annotation(Dialog(group = "Flow Coefficient"));
      parameter SI.Area Av(
        fixed=if CvData == CvTypes.Av then true else false,
        start=m_flow_nominal/(sqrt(d_nominal*dp_nominal))*valveCharacteristic(
            stemPosition_nominal)) = 0 "Av (metric) flow coefficient" 
       annotation(Dialog(group = "Flow Coefficient",
                         enable = (CvData==Modelica_Fluid.Types.CvTypes.Av)));
      parameter Real Kv = 0 "Kv (metric) flow coefficient [m3/h]" 
      annotation(Dialog(group = "Flow Coefficient",
                        enable = (CvData==CvTypes.Kv)));
      parameter Real Cv = 0 "Cv (US) flow coefficient [USG/min]" 
      annotation(Dialog(group = "Flow Coefficient",
                        enable = (CvData==CvTypes.Cv)));
      parameter SI.Pressure dp_nominal "Nominal pressure drop" 
      annotation(Dialog(group="Nominal operating point"));
      parameter Medium.MassFlowRate m_flow_nominal "Nominal mass flowrate" 
      annotation(Dialog(group="Nominal operating point"));
      parameter Medium.Density d_nominal=Medium.density_pTX(Medium.p_default, Medium.T_default, Medium.X_default)
        "Nominal inlet density" 
      annotation(Dialog(group="Nominal operating point"));
      parameter Real stemPosition_nominal=1 "Nominal stem position" 
      annotation(Dialog(group="Nominal operating point"));
      parameter Boolean CheckValve=false "Reverse flow stopped";

      replaceable function valveCharacteristic = 
          Modelica_Fluid.ControlValves.BaseClasses.ValveCharacteristics.linear 
        constrainedby
        Modelica_Fluid.ControlValves.BaseClasses.ValveCharacteristics.baseFun
        "Inherent flow characteristic" 
        annotation(choicesAllMatching=true);

      parameter Real delta=0.01 "Regularisation factor" 
                                                      annotation(Dialog(tab="Advanced"));

      constant SI.Area Kv2Av = 27.7e-6 "Conversion factor";
      constant SI.Area Cv2Av = 24.0e-6 "Conversion factor";

      Modelica.Blocks.Interfaces.RealInput stemPosition(min=-1e-10, max=1)
        "Stem position in the range 0-1" 
                                       annotation (Placement(transformation(
            origin={0,90},
            extent={{-20,-20},{20,20}},
            rotation=270), iconTransformation(
            extent={{-20,-20},{20,20}},
            rotation=270,
            origin={0,80})));

      parameter Real minStemPosition(
        min=0,
        max=0.1) = 0 "Minimum stemPosition (leckage flow to improve numerics)" 
      annotation(Dialog(tab="Advanced"));
      Real modifiedStemPosition
        "Modified, actually used stemPosition, so that the valve is not completely closed to improve numerics";

      annotation (
        Icon(coordinateSystem(
            preserveAspectRatio=true,
            extent={{-100,-100},{100,100}},
            grid={2,2}), graphics={
            Text(
              extent={{-150,-78},{150,-118}},
              textString="%name",
              lineColor={0,0,255}),
            Line(points={{0,52},{0,0}}, color={0,0,0}),
            Polygon(
              points={{-100,50},{-100,-50},{0,0},{-100,50}},
              lineColor={0,0,0},
              fillColor={255,255,255},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{100,50},{0,0},{100,-50},{100,50}},
              lineColor={0,0,0},
              fillPattern=FillPattern.Solid,
              fillColor={255,255,255}),
            Rectangle(
              extent={{-20,60},{20,52}},
              lineColor={0,0,0},
              fillColor={0,0,0},
              fillPattern=FillPattern.Solid),
            Polygon(
              points={{20,-36},{60,-50},{20,-64},{20,-36}},
              lineColor={0,128,255},
              smooth=Smooth.None,
              fillColor={0,128,255},
              fillPattern=FillPattern.Solid),
            Line(
              points={{30,-50},{-60,-50}},
              color={0,128,255},
              smooth=Smooth.None)}),
        Diagram(coordinateSystem(
            preserveAspectRatio=false,
            extent={{-100,-100},{100,100}},
            grid={2,2}), graphics),
        Documentation(info="<HTML>
<p>This is the base model for the <tt>ValveIncompressible</tt>, <tt>ValveVaporizing</tt>, and <tt>ValveCompressible</tt> valve models. The model is based on the IEC 534 / ISA S.75 standards for valve sizing.
<p>The model optionally supports reverse flow conditions (assuming symmetrical behaviour) or check valve operation, and has been suitably modified to avoid numerical singularities at zero pressure drop. 
<p><b>Modelling options</b></p>
<p>The following options are available to specify the valve flow coefficient in fully open conditions:
<ul><li><tt>CvData = Modelica_Fluid.Types.CvTypes.Av</tt>: the flow coefficient is given by the metric <tt>Av</tt> coefficient (m^2).
<li><tt>CvData = Modelica_Fluid.Types.CvTypes.Kv</tt>: the flow coefficient is given by the metric <tt>Kv</tt> coefficient (m^3/h).
<li><tt>CvData = Modelica_Fluid.Types.CvTypes.Cv</tt>: the flow coefficient is given by the US <tt>Cv</tt> coefficient (USG/min).
<li><tt>CvData = Modelica_Fluid.Types.CvTypes.OpPoint</tt>: the flow is computed from the nominal operating point specified by <tt>p_nominal</tt>, <tt>dp_nominal</tt>, <tt>m_flow_nominal</tt>, <tt>d_nominal</tt>, <tt>stemPosition_nominal</tt>.
</ul>
<p>The nominal pressure drop <tt>dp_nominal</tt> must always be specified; to avoid numerical singularities, the flow characteristic is modified for pressure drops less than <tt>b*dp_nominal</tt> (the default value is 1% of the nominal pressure drop). Increase this parameter if numerical problems occur in valves with very low pressure drops.
<p>If <tt>CheckValve</tt> is true, then the flow is stopped when the outlet pressure is higher than the inlet pressure; otherwise, reverse flow takes place. Use this option only when neede, as it increases the numerical complexity of the problem.
<p>The valve opening characteristic <tt>valveCharacteristic</tt>, linear by default, can be replaced by any user-defined function. Quadratic and equal percentage with customizable rangeability are already provided by the library.
</p>
 
<p>
The treatment of parameters <b>Kv</b> and <b>Cv</b> is
explained in detail in the 
<a href=\"Modelica://Modelica_Fluid.UsersGuide.ComponentDefinition.ValveCharacteristics\">Users Guide</a>.
</p>
 
</HTML>", revisions="<html>
<ul>
<li><i>2 Nov 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted from the ThermoPower library.</li>
</ul>
</html>"));
    protected
      function sqrtR = Utilities.regRoot(delta = delta*dp_nominal);

    initial equation
      if CvData == CvTypes.Kv then
        Av = Kv*Kv2Av "Unit conversion";
      elseif CvData == CvTypes.Cv then
        Av = Cv*Cv2Av "Unit conversion";
      end if;

    equation
      modifiedStemPosition = noEvent(if stemPosition > minStemPosition then 
        stemPosition else minStemPosition);
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
