package Turbomachinery 
  model Pump "Centrifugal pump with ideally controlled speed" 
    extends Interfaces.PartialPump;
    import Modelica.SIunits.Conversions.NonSIunits.*;
    parameter AngularVelocity_rpm N_const = N_nom "Constant rotational speed";
    Modelica.Blocks.Interfaces.RealInput N_in(redeclare type SignalType = 
          Modelica.SIunits.Conversions.NonSIunits.AngularVelocity_rpm) 
      "Prescribed rotational speed" 
      annotation (extent=[-36,34; -16,54],   rotation=-90);
  equation 
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
</HTML>",
        revisions="<html>
<ul>
<li><i>31 Oct 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
  end Pump;
  
  model PumpShaft "Centrifugal pump with mechanical connector for the shaft" 
    extends Interfaces.PartialPump;
    SI.Angle phi "Shaft angle";
    SI.AngularVelocity omega "Shaft angular velocity";
    Modelica.Mechanics.Rotational.Interfaces.Flange_a shaft 
    annotation (extent=[80,4; 110,32]);
  equation 
    phi = shaft.phi;
    omega = der(phi);
    N = Modelica.SIunits.Conversions.to_rpm(omega);
    W_single = omega*shaft.tau;
  annotation (
    Icon(Rectangle(extent=[60,26; 84,12], style(
            color=10,
            rgbcolor={95,95,95},
            gradient=2,
            fillColor=10,
            rgbfillColor={95,95,95}))),
    Diagram,
    Documentation(info="<HTML>
<p>This model describes a centrifugal pump (or a group of <tt>Np</tt> pumps in parallel) with a mechanical rotational connector for the shaft, to be used when the pump drive has to be modelled explicitly. In the case of <tt>Np</tt> pumps in parallel, the mechanical connector is relative to a single pump.
<p>The model extends <tt>PartialPump</tt>
 </HTML>",
       revisions="<html>
<ul>
<li><i>31 Oct 2005</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Model added to the Fluid library</li>
</ul>
</html>"));
  end PumpShaft;
  
end Turbomachinery;
