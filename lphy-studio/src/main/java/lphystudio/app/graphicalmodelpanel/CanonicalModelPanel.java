package lphystudio.app.graphicalmodelpanel;

import lphy.core.codebuilder.CanonicalCodeBuilder;
import lphy.core.logger.LoggerUtils;
import lphystudio.core.codecolorizer.DataModelCodeColorizer;

import javax.swing.*;
import javax.swing.text.BadLocationException;
import javax.swing.text.rtf.RTFEditorKit;
import java.awt.*;

/**
 * The panel to print the canonical model generated by the {@link GraphicalModelParserDictionary}.
 */
public class CanonicalModelPanel extends JComponent {
    GraphicalModelParserDictionary parserDictionary;
    JTextPane pane = new JTextPane();
    JScrollPane scrollPane;

    CanonicalCodeBuilder codeBuilder = new CanonicalCodeBuilder();

    public CanonicalModelPanel(GraphicalModelParserDictionary parserDictionary) {
        this.parserDictionary = parserDictionary;

        pane.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 12));
        pane.setBorder(BorderFactory.createEmptyBorder(10,10,10,10));
        pane.setEditable(false);
        pane.setEditorKit(new RTFEditorKit());

        JPanel noWrapPanel = new JPanel(new BorderLayout());
        noWrapPanel.add(pane);
        scrollPane = new JScrollPane(noWrapPanel);

        setText();

        BoxLayout boxLayout = new BoxLayout(this, BoxLayout.PAGE_AXIS);
        setLayout(boxLayout);
        add(scrollPane);

        parserDictionary.addGraphicalModelChangeListener(this::setText);
    }

    private void setText() {

        try {
            pane.getDocument().remove(0, pane.getDocument().getLength());
        } catch (BadLocationException e) {
            e.printStackTrace();
        }

        String text = codeBuilder.getCode(parserDictionary);

        System.out.println(text);

        if (text.length() > 0) {
            try {
                pane.setEditorKit(new RTFEditorKit());
                DataModelCodeColorizer codeColorizer = new DataModelCodeColorizer(parserDictionary, pane);
                codeColorizer.parse(text);

            }  catch (Exception e) {
                pane.setText(text);
                LoggerUtils.log.severe("CanonicalModelPanel failed with exception: " + e.getMessage());
            }
        }
    }
}
