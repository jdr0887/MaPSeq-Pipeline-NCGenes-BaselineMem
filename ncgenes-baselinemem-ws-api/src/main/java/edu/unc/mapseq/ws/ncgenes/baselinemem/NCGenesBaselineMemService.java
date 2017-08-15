package edu.unc.mapseq.ws.ncgenes.baseline;

import javax.jws.WebMethod;
import javax.jws.WebParam;
import javax.jws.WebService;
import javax.jws.soap.SOAPBinding;
import javax.jws.soap.SOAPBinding.Style;
import javax.jws.soap.SOAPBinding.Use;
import javax.ws.rs.Consumes;
import javax.ws.rs.GET;
import javax.ws.rs.Path;
import javax.ws.rs.PathParam;
import javax.ws.rs.Produces;
import javax.ws.rs.core.MediaType;
import javax.xml.ws.BindingType;

import org.renci.vcf.VCFResult;

@BindingType(value = javax.xml.ws.soap.SOAPBinding.SOAP11HTTP_BINDING)
@WebService(targetNamespace = "http://baseline.ncgenes.ws.mapseq.unc.edu", serviceName = "NCGenesBaselineMemService", portName = "NCGenesBaselineMemPort")
@SOAPBinding(style = Style.DOCUMENT, use = Use.LITERAL, parameterStyle = SOAPBinding.ParameterStyle.WRAPPED)
@Path("/NCGenesBaselineMemService/")
@Consumes(MediaType.APPLICATION_JSON)
@Produces(MediaType.APPLICATION_JSON)
public interface NCGenesBaselineMemService {

    @GET
    @Path("/lookupQuantificationResults/{sampleId}")
    @WebMethod
    public QualityControlInfo lookupQuantificationResults(@PathParam("sampleId") @WebParam(name = "sampleId") Long sampleId);

    @GET
    @Path("/lookupIdentityInfoFromVCF/{sampleId}")
    @WebMethod
    public VCFResult lookupIdentityInfoFromVCF(@PathParam("sampleId") @WebParam(name = "sampleId") Long sampleId);

}
