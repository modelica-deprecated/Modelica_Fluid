package Pipes "Lumped, distributed and thermal pipe components" 
  extends Icons.VariantLibrary;
  
  model AsymmetricDistributedPipe 
    "Distributed pipe model using asymmetric discretizations" 
    // Taken one-to-one from Modelica_Fluid
    
    // Covered herein:
    //   Friction pressure losses at ports
    
    extends Interfaces.PartialComponent;
    
    extends FluidInterface.PartialAsymmetricDistributedPipe;
    
    annotation (Documentation(info="<html>
<p>Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow\">PartialDistributedFlow</a>. Source terms in the mass balances are set to zero. The total volume is a parameter. The additional component <tt>heatTransfer</tt> specifies the source term <tt>Qs_flow</tt> in the energy balance. The default component uses a constant coefficient of heat transfer to model convective heat transfer between segment boundary (<tt>thermalPort</tt>) and the bulk flow. The <tt>heatTransfer</tt> model is replaceable and can be exchanged with any model extended from <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer\">PartialPipeHeatTransfer</a>.</p>
<p>Pressure drop correlations (algebraic and possibly non-linear flow model) correlate the pressure in the first control volume with the pressure in port_a and the pressures of port_b and the nth control volume, respectively.</p>
 
</html>",   revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"),
      Diagram,
      Icon(Rectangle(extent=[-100,40; 100,-40], style(
            color=69,
            gradient=2,
            fillColor=69)),
        Rectangle(extent=[-72,40; -24,-40], style(color=0, rgbcolor={0,0,0})),
        Rectangle(extent=[-24,40; 24,-40], style(color=0, rgbcolor={0,0,0})),
        Rectangle(extent=[24,40; 72,-40], style(color=0, rgbcolor={0,0,0})),
        Rectangle(extent=[72,40; 120,-40], style(color=0, rgbcolor={0,0,0}))));
  equation 
    assert(FluidDiscretization.isSymmetric == false,
      "This model uses a pipe component with Port_A Port_B but the chosen discretization is symmetric. Please use a Port_B Port_B implementation.");
    
  end AsymmetricDistributedPipe;
  
  model SymmetricDistributedPipe 
    "Distributed pipe model using symmetric discretizations" 
    // Taken one-to-one from Modelica_Fluid
    
    // Covered herein:
    //   Friction pressure losses at ports
    
    extends Interfaces.PartialComponent;
    
    extends FluidInterface.PartialSymmetricDistributedPipe;
    
    annotation (Documentation(info="<html>
<p>Distributed pipe model based on <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.PartialDistributedFlow\">PartialDistributedFlow</a>. Source terms in the mass balances are set to zero. The total volume is a parameter. The additional component <tt>heatTransfer</tt> specifies the source term <tt>Qs_flow</tt> in the energy balance. The default component uses a constant coefficient of heat transfer to model convective heat transfer between segment boundary (<tt>thermalPort</tt>) and the bulk flow. The <tt>heatTransfer</tt> model is replaceable and can be exchanged with any model extended from <a href=\"Modelica:Modelica_Fluid.Pipes.BaseClasses.HeatTransfer.PartialPipeHeatTransfer\">PartialPipeHeatTransfer</a>.</p>
<p>Pressure drop correlations (algebraic and possibly non-linear flow model) correlate the pressure in the first control volume with the pressure in port_a and the pressures of port_b and the nth control volume, respectively.</p>
 
</html>",   revisions="<html>
<ul>
<li><i>04 Mar 2006</i>
    by Katrin Pr&ouml;l&szlig;:<br>
       Model added to the Fluid library</li>
</ul>
</html>"),
      Diagram,
      Icon(Rectangle(extent=[-100,40; 100,-40], style(
            color=69,
            gradient=2,
            fillColor=69)),
        Rectangle(extent=[-96,40; -72,-40], style(color=0, rgbcolor={0,0,0})),
        Rectangle(extent=[-72,40; -24,-40], style(color=0, rgbcolor={0,0,0})),
        Rectangle(extent=[24,40; 72,-40], style(color=0, rgbcolor={0,0,0})),
        Rectangle(extent=[-24,40; 24,-40], style(color=0, rgbcolor={0,0,0})),
        Rectangle(extent=[72,40; 96,-40], style(color=0, rgbcolor={0,0,0}))));
  equation 
    assert(FluidDiscretization.isSymmetric,
      "This model uses a pipe component with Port_B Port_B but the chosen discretization is asymmetric. Please use a Port_A Port_B implementation.");
    
  end SymmetricDistributedPipe;
end Pipes;
