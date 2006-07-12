package HeatTransfer 
  model PipeHT_constAlpha 
    extends Modelica_Fluid.Interfaces.HeatTransfer.PartialPipeHeatTransfer;
    parameter SI.CoefficientOfHeatTransfer alpha0=200;
    annotation(structurallyIncomplete=true, Documentation(info="<html>
Simple heat transfer correlation with constant heat transfer coefficient, used as default component in <a href=\"Modelica:Modelica_Fluid.Components.Pipes.DistributedPipe_thermal\">DistributedPipe_thermal</a>.
</html>"));
  equation 
    for i in 1:n loop
      thermalPort[i].Q_flow=alpha0*A_h/n*(thermalPort[i].T-T[i]);
    end for;
    thermalPort.Q_flow=Q_flow;
  end PipeHT_constAlpha;
end HeatTransfer;
