package browser.files;

import browser.Cache;
import browser.ICONS_MANAGER;
import browser.LABELS;
import browser.imageitems.listitems.FileItem;
import browser.imageitems.listitems.FileItemFolder;
import browser.imageitems.listitems.FileImageItem;
import browser.imageitems.SelFileItem;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.regex.Pattern;
import javax.swing.AbstractListModel;
import javax.swing.JLabel;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.text.BadLocationException;
import javax.swing.text.Document;

/**
 * Manages filtering of list model
 */
public class FilterFilesModel extends AbstractListModel implements DocumentListener {

    protected final static FileItem PARENT_DIRECTORY = createSuitableFileItem(new File(".." + System.getProperty("file.separator")));
    protected final static String PATTERNS[] = {"*", "?"};
    protected final static String STATIC_PATTERNS[] = {"."};
    protected JLabel filteringLabel;
    protected FileBrowser fileBrowser;
    protected List<FileItem> list;
    protected List<FileItem> filteredList;
    protected String lastFilter = "";
    protected static Cache cache = new Cache(ICONS_MANAGER.CACHE_SIZE_BROWSER);

    public FilterFilesModel(String work_dir) {
        super();

        list = new ArrayList<FileItem>();
        filteredList = new ArrayList<FileItem>();

        fileBrowser = new FileBrowser(work_dir);

        refreshList();
    }

    public void setFilteringLabel(JLabel filteringLabel) {
        this.filteringLabel = filteringLabel;
    }

    public void addElement(Object element) {
        list.add((FileItem) element);
    }

    public int getTotalSize() {
        return list.size();
    }

    public int getSize() {
        return filteredList.size();
    }

    public Object getElementAt(int index) {
        Object returnValue = null;

        if (index < filteredList.size()) {
            returnValue = filteredList.get(index);
        }

        return returnValue;
    }

    private static String fixRegularExpression(String regexp, String patterns[], String prefix) {
        StringBuffer fixed = new StringBuffer(regexp);

        for (int i = 0; i < patterns.length; i++) {
            int index = 0;
            while ((index = fixed.indexOf(patterns[i], index)) >= 0) {
                fixed.insert(index, prefix);
                index += 1 + patterns[i].length();
            }
        }

        return fixed.toString();
    }

    private static String fixRegularExpression(String regexp) {
        String fixed = fixRegularExpression(regexp, STATIC_PATTERNS, "\\");
        fixed = fixRegularExpression(fixed, PATTERNS, ".");

        return fixed;
    }

    // DocumentListener Methods
    public void insertUpdate(DocumentEvent event) {
        Document doc = event.getDocument();
        try {
            lastFilter = doc.getText(0, doc.getLength());
            filterList();
        } catch (BadLocationException ble) {
            System.err.println("Bad location: " + ble);
        }
    }

    public void removeUpdate(DocumentEvent event) {
        Document doc = event.getDocument();
        try {
            lastFilter = doc.getText(0, doc.getLength());
            filterList();
        } catch (BadLocationException ble) {
            System.err.println("Bad location: " + ble);
        }
    }

    public void changedUpdate(DocumentEvent event) {
    }

    public String getCurrentDirectory() {
        return fileBrowser.getCurrentDirectory();
    }

    public File getCurrentRoot() {
        return fileBrowser.getCurrentRoot();
    }

    public void changeDirectory(int index) {
        if (index > 0) {
            fileBrowser.changeDirectory(((FileItem) getElementAt(index)).getFile());
        } else {
            goParent();
        }

        refreshList();
    }

    public void changeDirectory(File newDir) {
        fileBrowser.changeDirectory(newDir);
        refreshList();
    }

    public void goParent() {
        fileBrowser.goParent();
        refreshList();
    }

    public void refreshCurrentDirectory() {
        fileBrowser.refresh();
        refreshList();
    }

    protected void refreshList() {
        buildList();
        fireContentsChanged(this, 0, getSize());
    }

    protected void buildList() {
        list.clear();   // Remove all.

        for (File file : fileBrowser.getFiles()) {
            addElement(createSuitableFileItem(file));
        }

        filterList();
    }

    protected static FileItem createSuitableFileItem(File file) {
        FileItem fileItem;

        if (file.isDirectory()) {
            fileItem = new FileItemFolder(file);
        } else if (FileBrowser.isBasicImageType(file) || FileBrowser.isSpiderImage(file)) {
            fileItem = new FileImageItem(file, cache);
        } else if (FileBrowser.isSelFile(file)) {
            fileItem = new SelFileItem(file, cache);
        } else {
            fileItem = new FileItem(file);
        }

        return fileItem;
    }

    private void filterList() {
        filteredList.clear();   // Clear list.
        filteredList.add(PARENT_DIRECTORY); // ".." directory will be always present.

        // Gets all tokens from filter.
        StringTokenizer st = new StringTokenizer(lastFilter, " ");
        String filters[] = new String[st.countTokens()];
        for (int i = 0; i < filters.length; i++) {
            filters[i] = st.nextToken();
        }

        for (FileItem file : list) {
            if (filters.length == 0 || filter(file.getFile().getName(), filters)) {
                filteredList.add(file);
            }
        }

        fireContentsChanged(this, 0, getSize());

        // Sets label to show filtering message (if any).
        if (filteringLabel != null) {
            int showingItems = getSize() - 1;
            int totalItems = getTotalSize();

            filteringLabel.setIcon(showingItems != totalItems ? ICONS_MANAGER.FILTERING_ALERT : null);
            filteringLabel.setText(LABELS.LABEL_FILES_SHOWING(showingItems, totalItems));
            filteringLabel.updateUI();
        }
    }

    private static boolean filter(String file, String filters[]) {
        for (int i = 0; i < filters.length; i++) {
            if (filter(file, filters[i])) {
                return true;
            }
        }

        return false;
    }

    private static boolean filter(String file, String search) {
        String fixedExpresion = fixRegularExpression(search);

        return Pattern.compile(fixedExpresion, Pattern.CASE_INSENSITIVE).matcher(file).matches();
    }
}
