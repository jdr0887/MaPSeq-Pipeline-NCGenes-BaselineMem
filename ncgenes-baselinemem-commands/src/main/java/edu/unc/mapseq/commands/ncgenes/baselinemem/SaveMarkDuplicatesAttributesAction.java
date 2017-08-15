package edu.unc.mapseq.commands.ncgenes.baselinemem;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.karaf.shell.api.action.Action;
import org.apache.karaf.shell.api.action.Command;
import org.apache.karaf.shell.api.action.Option;
import org.apache.karaf.shell.api.action.lifecycle.Reference;
import org.apache.karaf.shell.api.action.lifecycle.Service;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import edu.unc.mapseq.commons.ncgenes.baselinemem.SaveMarkDuplicatesAttributesRunnable;
import edu.unc.mapseq.dao.MaPSeqDAOBeanService;
import edu.unc.mapseq.dao.model.WorkflowRunAttempt;

@Command(scope = "ncgenes-baselinemem", name = "save-mark-duplicates-attributes", description = "Save MarkDuplicates Attributes")
@Service
public class SaveMarkDuplicatesAttributesAction implements Action {

    private static final Logger logger = LoggerFactory.getLogger(SaveMarkDuplicatesAttributesAction.class);

    @Reference
    private MaPSeqDAOBeanService maPSeqDAOBeanService;

    @Option(name = "--workflowRunAttemptId", description = "WorkflowRunAttempt Identifier", required = true, multiValued = false)
    private Long workflowRunAttemptId;

    @Override
    public Object execute() throws Exception {
        logger.debug("ENTERING execute()");
        ExecutorService es = Executors.newSingleThreadExecutor();
        WorkflowRunAttempt attempt = maPSeqDAOBeanService.getWorkflowRunAttemptDAO().findById(workflowRunAttemptId);
        SaveMarkDuplicatesAttributesRunnable runnable = new SaveMarkDuplicatesAttributesRunnable(maPSeqDAOBeanService, attempt);
        es.submit(runnable);
        es.shutdown();
        return null;
    }

    public Long getWorkflowRunAttemptId() {
        return workflowRunAttemptId;
    }

    public void setWorkflowRunAttemptId(Long workflowRunAttemptId) {
        this.workflowRunAttemptId = workflowRunAttemptId;
    }

}
