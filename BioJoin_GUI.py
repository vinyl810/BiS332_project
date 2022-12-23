# import requirements

from PySide6.QtCore import (
    QDate,
    QMetaObject,
    QRect,
    Qt,
    Slot,
)

from PySide6.QtGui import (
    QAction,
    QIntValidator,
    QValidator,
)

from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QComboBox,
    QRadioButton,
    QDateEdit,
    QDialog,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMenu,
    QMenuBar,
    QPlainTextEdit,
    QPushButton,
    QStackedWidget,
    QStatusBar,
    QTableWidget,
    QTableWidgetItem,
    QTextBrowser,
    QVBoxLayout,
    QWidget,
    QMessageBox,
)

from pandas import DataFrame

from query import sendQuery
from geneVis import display_snp
import genetic_distance as GD

import os.path
from math import isnan

HOME = """\
==============================================================
.______    __    ______          __    ______    __  .__   __. 
|   _  \  |  |  /  __  \        |  |  /  __  \  |  | |  \ |  | 
|  |_)  | |  | |  |  |  |       |  | |  |  |  | |  | |   \|  | 
|   _  <  |  | |  |  |  | .--.  |  | |  |  |  | |  | |  . `  | 
|  |_)  | |  | |  `--'  | |  `--'  | |  `--'  | |  | |  |\   | 
|______/  |__|  \______/   \______/   \______/  |__| |__| \__| 

            Made by | 18 Daewon K. 20 Haejoon P. 21 Gyuchan A.
==============================================================
"""
INSERT = """\
 __  .__   __.      _______. _______ .______     .___________.
|  | |  \ |  |     /       ||   ____||   _  \    |           |
|  | |   \|  |    |   (----`|  |__   |  |_)  |   `---|  |----`
|  | |  . `  |     \   \    |   __|  |      /        |  |     
|  | |  |\   | .----)   |   |  |____ |  |\  \----.   |  |     
|__| |__| \__| |_______/    |_______|| _| `._____|   |__|     
"""
UPDATE = """\
 __    __  .______    _______       ___   .___________. _______ 
|  |  |  | |   _  \  |       \     /   \  |           ||   ____|
|  |  |  | |  |_)  | |  .--.  |   /  ^  \ `---|  |----`|  |__   
|  |  |  | |   ___/  |  |  |  |  /  /_\  \    |  |     |   __|  
|  `--'  | |  |      |  '--'  | /  _____  \   |  |     |  |____ 
 \______/  | _|      |_______/ /__/     \__\  |__|     |_______|
"""
DELETE = """\
 _______   _______  __       _______ .___________. _______ 
|       \ |   ____||  |     |   ____||           ||   ____|
|  .--.  ||  |__   |  |     |  |__   `---|  |----`|  |__   
|  |  |  ||   __|  |  |     |   __|      |  |     |   __|  
|  '--'  ||  |____ |  `----.|  |____     |  |     |  |____ 
|_______/ |_______||_______||_______|    |__|     |_______|
"""
SEARCH = """\
     _______. _______     ___      .______        ______  __    __  
    /       ||   ____|   /   \     |   _  \      /      ||  |  |  | 
   |   (----`|  |__     /  ^  \    |  |_)  |    |  ,----'|  |__|  | 
    \   \    |   __|   /  /_\  \   |      /     |  |     |   __   | 
.----)   |   |  |____ /  _____  \  |  |\  \----.|  `----.|  |  |  | 
|_______/    |_______/__/     \__\ | _| `._____| \______||__|  |__| 
"""
DISTANCE = """\
 _______   __       _______.___________.    ___      .__   __.   ______  _______ 
|       \ |  |     /       |           |   /   \     |  \ |  |  /      ||   ____|
|  .--.  ||  |    |   (----`---|  |----`  /  ^  \    |   \|  | |  ,----'|  |__   
|  |  |  ||  |     \   \       |  |      /  /_\  \   |  . `  | |  |     |   __|  
|  '--'  ||  | .----)   |      |  |     /  _____  \  |  |\   | |  `----.|  |____ 
|_______/ |__| |_______/       |__|    /__/     \__\ |__| \__|  \______||_______|
"""

INT_MAX = 2147483647
TIMEOUT_MS = 5000

# TODO: bad OOP
# TODO: too many magic numbers
# TODO: meta-approach for simple code (getattr, setattr)
FIELD_SNP = ("snp_id", "chromo", "pos", "neighbor_gene", "anc_allele", "minor_allele")
LABEL_SNP = (
    "SNP ID",
    "Chromosome",
    "Position",
    "Neighbor gene",
    "Ancestor allele",
    "Minor allele",
)
PKIDX_SNP = 0
FIELD_GENE = (
    "tax_id",
    "gene_id",
    "symbol",
    "synonym",
    "chromo",
    "map_loc",
    "description",
    "type_gene",
    "mod_date",
)
LABEL_GENE = (
    "Tax ID",
    "Gene ID",
    "Symbol",
    "Synonym",
    "Chromosome",
    "Map location",
    "Description",
    "Gene type",
    "Modification date",
)
PKIDX_GENE = 1
FIELD_DISEASE = ("disease_id", "name")
LABEL_DISEASE = ("Disease ID", "Name")
PKIDX_DISEASE = 0

opened_windows = set()
log_window = None
log_list = list()


def getResponse(query: str):
    log_list.append("[query]\n" + query)
    result = sendQuery(query)
    log_list.append("[reply]\n" + str(result))
    if result is None:
        return "Successfully Done."
    else:
        return result


def exportLog():
    with open("BioJoin.log", "w") as logfile:
        logfile.write(logToStr())


def logToStr():
    global log_list
    return "\n\n".join(log_list)


def fieldCheckValid(field: QWidget):
    if hasattr(field, "text"):
        text = field.text()
    elif hasattr(field, "currentText"):
        text = field.currentText()
    else:
        raise TypeError("Unavailable type:", type(field).__name__())

    if hasattr(field, "validator"):
        if (
            field.validator() is None
            or field.validator().validate(text, 0)[0] == QValidator.Acceptable
        ):
            return True
        else:
            return False
    elif hasattr(field, "validate"):
        if field.validate(text, 0)[0] == QValidator.Acceptable:
            return True
        else:
            return False
    else:
        return True


def formToTuple(form: QFormLayout, allowEmptyForInt=False):
    text_list = list()
    for i in range(form.rowCount()):
        field = form.itemAt(i, QFormLayout.FieldRole).widget()
        if hasattr(field, "text"):
            text = field.text()
        elif hasattr(field, "currentText"):
            text = field.currentText()
        else:
            raise TypeError("Unavailable type:", type(field).__name__())

        if (
            hasattr(field, "validator")
            and not field.validator() is None
            and isinstance(field.validator(), QIntValidator)
        ):
            try:
                text_list.append(int(text))
            except ValueError as e:
                if allowEmptyForInt:
                    text_list.append("")
                else:
                    raise e
        else:
            text_list.append(text)

    return tuple(text_list)


def searchResultToForm(_form: QFormLayout, searchResult: DataFrame):
    assert _form.rowCount() == len(searchResult.columns)
    for (i, value) in enumerate(searchResult.values.flatten().tolist()):
        field = _form.itemAt(i, QFormLayout.FieldRole).widget()
        if isinstance(field, QLineEdit):
            field.setText(str(value))
        elif isinstance(field, QComboBox):
            field.setCurrentText(value)
        elif isinstance(field, QDateEdit):
            field.setDate(QDate(value))
        else:
            raise TypeError("Unavailable type:", type(field).__name__())


def searchResultToTable(_table: QTableWidget, searchResult: DataFrame, labels):
    (rowCount, colCount) = searchResult.shape
    assert rowCount > 0
    assert len(labels) == colCount
    _table.clear()
    _table.setRowCount(rowCount)
    _table.setColumnCount(colCount)
    _table.setHorizontalHeaderLabels(labels)
    _table.scrollToTop()
    for i in range(rowCount):
        for j in range(colCount):
            data = searchResult.iat[i, j]
            if isinstance(data, float):
                if isnan(data):
                    data = None
                else:
                    data = int(data)
            _table.setItem(i, j, QTableWidgetItem(str(data)))


def quoteWrappedStr(value):
    if isinstance(value, str):
        return "".join(["'", value, "'"])
    else:
        return str(value)


class AlphanumericValidator(QValidator):
    def validate(self, input: str, pos: int):
        if input.replace("-", "").isalnum():
            return (QValidator.Acceptable, input, pos)
        else:
            if len(input) == 0:
                return (QValidator.Intermediate, input, pos)
            return (QValidator.Invalid, input, pos)

    def fixup(self, input: str):
        pass


class ChromosomeValidator(QValidator):
    def validate(self, input: str, pos: int):
        try:
            num = int(input)
            if 1 <= num <= 22:
                return (QValidator.Acceptable, input, pos)
            if num <= 29:
                return (QValidator.Intermediate, input, pos)
            return (QValidator.Invalid, input, pos)
        except ValueError:
            if input.upper() == "X" or input.upper() == "Y":
                return (QValidator.Acceptable, input.upper(), pos)
            if len(input) == 0:
                return (QValidator.Intermediate, input, pos)
            return (QValidator.Invalid, input, pos)

    def fixup(self, input: str):
        pass


class InsertSnpWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUI()

    def setupUI(self):
        self.form = QFormLayout(self)
        self.form.setContentsMargins(0, 0, 0, 0)

        self.label_insert_snp_id = QLabel(self)
        self.form.setWidget(0, QFormLayout.LabelRole, self.label_insert_snp_id)
        self.label_insert_snp_id.setText("SNP ID")

        self.label_insert_snp_chromo = QLabel(self)
        self.form.setWidget(1, QFormLayout.LabelRole, self.label_insert_snp_chromo)
        self.label_insert_snp_chromo.setText("Chromosome")

        self.label_insert_snp_pos = QLabel(self)
        self.form.setWidget(2, QFormLayout.LabelRole, self.label_insert_snp_pos)
        self.label_insert_snp_pos.setText("Position")

        self.label_insert_snp_neighbor = QLabel(self)
        self.form.setWidget(3, QFormLayout.LabelRole, self.label_insert_snp_neighbor)
        # This is the longest label.
        self.label_insert_snp_neighbor.setText("Neighbor gene\t\t")

        self.label_insert_snp_anc = QLabel(self)
        self.form.setWidget(4, QFormLayout.LabelRole, self.label_insert_snp_anc)
        self.label_insert_snp_anc.setText("Ancestor allele")

        self.label_insert_snp_minor = QLabel(self)
        self.form.setWidget(5, QFormLayout.LabelRole, self.label_insert_snp_minor)
        self.label_insert_snp_minor.setText("Minor allele")

        self.lineEdit_insert_snp_id = QLineEdit(self)
        self.form.setWidget(0, QFormLayout.FieldRole, self.lineEdit_insert_snp_id)
        self.lineEdit_insert_snp_id.setValidator(QIntValidator(1, INT_MAX))
        self.lineEdit_insert_snp_id.setPlaceholderText("positive integer")

        self.lineEdit_insert_snp_chromo = QLineEdit(self)
        self.form.setWidget(1, QFormLayout.FieldRole, self.lineEdit_insert_snp_chromo)
        self.lineEdit_insert_snp_chromo.setValidator(ChromosomeValidator())
        self.lineEdit_insert_snp_chromo.setPlaceholderText("1~22 or X/Y")

        self.lineEdit_insert_snp_pos = QLineEdit(self)
        self.form.setWidget(2, QFormLayout.FieldRole, self.lineEdit_insert_snp_pos)
        self.lineEdit_insert_snp_pos.setValidator(QIntValidator(1, INT_MAX))
        self.lineEdit_insert_snp_pos.setPlaceholderText("positive integer")

        self.lineEdit_insert_snp_neighbor = QLineEdit(self)
        self.form.setWidget(3, QFormLayout.FieldRole, self.lineEdit_insert_snp_neighbor)
        self.lineEdit_insert_snp_neighbor.setValidator(AlphanumericValidator())
        self.lineEdit_insert_snp_neighbor.setPlaceholderText("alphanumeric and hyphen")

        self.comboBox_insert_snp_anc = QComboBox(self)
        self.form.setWidget(4, QFormLayout.FieldRole, self.comboBox_insert_snp_anc)
        self.comboBox_insert_snp_anc.addItem("A")
        self.comboBox_insert_snp_anc.addItem("C")
        self.comboBox_insert_snp_anc.addItem("G")
        self.comboBox_insert_snp_anc.addItem("T")

        self.comboBox_insert_snp_minor = QComboBox(self)
        self.form.setWidget(5, QFormLayout.FieldRole, self.comboBox_insert_snp_minor)
        self.comboBox_insert_snp_minor.addItem("A")
        self.comboBox_insert_snp_minor.addItem("C")
        self.comboBox_insert_snp_minor.addItem("G")
        self.comboBox_insert_snp_minor.addItem("T")


class InsertGeneWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUI()

    def setupUI(self):
        self.form = QFormLayout(self)
        self.form.setContentsMargins(0, 0, 0, 0)

        self.label_insert_gene_tax_id = QLabel(self)
        self.form.setWidget(0, QFormLayout.LabelRole, self.label_insert_gene_tax_id)
        self.label_insert_gene_tax_id.setText("Tax ID")

        self.label_insert_gene_id = QLabel(self)
        self.form.setWidget(1, QFormLayout.LabelRole, self.label_insert_gene_id)
        self.label_insert_gene_id.setText("Gene ID")

        self.label_insert_gene_symbol = QLabel(self)
        self.form.setWidget(2, QFormLayout.LabelRole, self.label_insert_gene_symbol)
        self.label_insert_gene_symbol.setText("Symbol")

        self.label_insert_gene_synonym = QLabel(self)
        self.form.setWidget(3, QFormLayout.LabelRole, self.label_insert_gene_synonym)
        self.label_insert_gene_synonym.setText("Synonym")

        self.label_insert_gene_chromo = QLabel(self)
        self.form.setWidget(4, QFormLayout.LabelRole, self.label_insert_gene_chromo)
        self.label_insert_gene_chromo.setText("Chromosome")

        self.label_insert_gene_map_loc = QLabel(self)
        self.form.setWidget(5, QFormLayout.LabelRole, self.label_insert_gene_map_loc)
        # This is the second-longest label. (Need for align in search_gene)
        self.label_insert_gene_map_loc.setText("Map location\t\t")

        self.label_insert_gene_desc = QLabel(self)
        self.form.setWidget(6, QFormLayout.LabelRole, self.label_insert_gene_desc)
        self.label_insert_gene_desc.setText("Description")

        self.label_insert_gene_type = QLabel(self)
        self.form.setWidget(7, QFormLayout.LabelRole, self.label_insert_gene_type)
        self.label_insert_gene_type.setText("Gene type")

        self.label_insert_gene_date = QLabel(self)
        self.form.setWidget(8, QFormLayout.LabelRole, self.label_insert_gene_date)
        # This is the longest label.
        self.label_insert_gene_date.setText("Modification date\t")

        self.lineEdit_insert_gene_tax_id = QLineEdit(self)
        self.form.setWidget(0, QFormLayout.FieldRole, self.lineEdit_insert_gene_tax_id)
        self.lineEdit_insert_gene_tax_id.setValidator(QIntValidator(1, INT_MAX))
        self.lineEdit_insert_gene_tax_id.setPlaceholderText("positive integer")

        self.lineEdit_insert_gene_id = QLineEdit(self)
        self.form.setWidget(1, QFormLayout.FieldRole, self.lineEdit_insert_gene_id)
        self.lineEdit_insert_gene_id.setValidator(QIntValidator(1, INT_MAX))
        self.lineEdit_insert_gene_id.setPlaceholderText("positive integer")

        self.lineEdit_insert_gene_symbol = QLineEdit(self)
        self.form.setWidget(2, QFormLayout.FieldRole, self.lineEdit_insert_gene_symbol)
        self.lineEdit_insert_gene_symbol.setValidator(AlphanumericValidator())
        self.lineEdit_insert_gene_symbol.setPlaceholderText("alphanumeric and hyphen")

        self.lineEdit_insert_gene_synonym = QLineEdit(self)
        self.form.setWidget(3, QFormLayout.FieldRole, self.lineEdit_insert_gene_synonym)
        # too long
        self.lineEdit_insert_gene_synonym.setPlaceholderText(
            "hover mouse to see format"
        )
        self.lineEdit_insert_gene_synonym.setToolTip(
            "Gene symbols separated by bar(|)\n(example) A1B2|C3D4|E5"
        )

        self.lineEdit_insert_gene_chromo = QLineEdit(self)
        self.form.setWidget(4, QFormLayout.FieldRole, self.lineEdit_insert_gene_chromo)
        self.lineEdit_insert_gene_chromo.setValidator(ChromosomeValidator())
        self.lineEdit_insert_gene_chromo.setPlaceholderText("1~22, X/Y")

        # TODO
        self.lineEdit_insert_gene_map_loc = QLineEdit(self)
        self.form.setWidget(5, QFormLayout.FieldRole, self.lineEdit_insert_gene_map_loc)
        self.lineEdit_insert_gene_map_loc.setPlaceholderText("??????")

        self.lineEdit_insert_gene_desc = QLineEdit(self)
        self.form.setWidget(6, QFormLayout.FieldRole, self.lineEdit_insert_gene_desc)
        self.lineEdit_insert_gene_desc.setPlaceholderText("description (string)")

        self.lineEdit_insert_gene_type = QLineEdit(self)
        self.form.setWidget(7, QFormLayout.FieldRole, self.lineEdit_insert_gene_type)
        self.lineEdit_insert_gene_type.setPlaceholderText("gene type (string)")

        self.dateEdit_insert_gene_date = QDateEdit(self)
        self.form.setWidget(8, QFormLayout.FieldRole, self.dateEdit_insert_gene_date)
        self.dateEdit_insert_gene_date.setDisplayFormat("yyyy-MM-dd")


class InsertDiseaseWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUI()

    def setupUI(self):
        self.form = QFormLayout(self)
        self.form.setContentsMargins(0, 0, 0, 0)

        self.label_insert_disease_id = QLabel(self)
        self.form.setWidget(0, QFormLayout.LabelRole, self.label_insert_disease_id)
        # This is the longest label.
        self.label_insert_disease_id.setText("Disease ID\t\t")

        self.label_insert_disease_name = QLabel(self)
        self.form.setWidget(1, QFormLayout.LabelRole, self.label_insert_disease_name)
        self.label_insert_disease_name.setText("Name")

        self.lineEdit_insert_disease_id = QLineEdit(self)
        self.form.setWidget(0, QFormLayout.FieldRole, self.lineEdit_insert_disease_id)
        self.lineEdit_insert_disease_id.setValidator(QIntValidator(1, INT_MAX))
        self.lineEdit_insert_disease_id.setPlaceholderText("positive integer")

        self.lineEdit_insert_disease_name = QLineEdit(self)
        self.form.setWidget(1, QFormLayout.FieldRole, self.lineEdit_insert_disease_name)
        self.lineEdit_insert_disease_name.setPlaceholderText("disease name (string)")


class MainWindow(QMainWindow):
    def __init__(self):
        global opened_windows
        super().__init__()
        self.dataBuffer = None
        self.setupUI()
        opened_windows.add(self)

    def setupUI(self):
        self.setFixedSize(800, 600)
        self.setWindowTitle("BioJoin")

        # region basic frames
        self.centralwidget = QWidget(self)
        self.setCentralWidget(self.centralwidget)

        self.stackedWidget = QStackedWidget(self.centralwidget)
        self.stackedWidget.setGeometry(QRect(0, 0, 800, 600))

        self.menubar = QMenuBar(self)
        self.setMenuBar(self.menubar)

        self.statusbar = QStatusBar(self)
        self.setStatusBar(self.statusbar)
        # endregion

        # region menubar
        self.menuMenu = QMenu(self.menubar)
        self.menubar.addAction(self.menuMenu.menuAction())
        self.menuMenu.setTitle("Menu")

        self.menuLog = QMenu(self.menubar)
        self.menubar.addAction(self.menuLog.menuAction())
        self.menuLog.setTitle("Log")

        self.menuAdvanced = QMenu(self.menubar)
        self.menubar.addAction(self.menuAdvanced.menuAction())
        self.menuAdvanced.setTitle("Advanced")
        # endregion

        # region menu
        self.actionHome = QAction(self)
        self.actionHome.setObjectName("actionHome")
        self.menuMenu.addAction(self.actionHome)
        self.actionHome.setShortcut("Ctrl+H")
        self.actionHome.setText("Home")

        self.actionOpen = QAction(self)
        self.actionOpen.setObjectName("actionOpen")
        self.menuMenu.addAction(self.actionOpen)
        self.actionOpen.setShortcut("Ctrl+N")
        self.actionOpen.setText("Open new window")

        self.actionClose = QAction(self)
        self.actionClose.setObjectName("actionClose")
        self.menuMenu.addAction(self.actionClose)
        self.actionClose.setShortcut("Esc")
        self.actionClose.setText("Close this window")

        self.actionShowLog = QAction(self)
        self.actionShowLog.setObjectName("actionShowLog")
        self.menuLog.addAction(self.actionShowLog)
        self.actionShowLog.setShortcut("Ctrl+L")
        self.actionShowLog.setText("Show log")

        self.actionExportLog = QAction(self)
        self.actionExportLog.setObjectName("actionExportLog")
        self.menuLog.addAction(self.actionExportLog)
        self.actionExportLog.setText("Export log")

        self.actionSQL = QAction(self)
        self.actionSQL.setObjectName("actionSQL")
        self.menuAdvanced.addAction(self.actionSQL)
        self.actionSQL.setText("SQL query...")
        # endregion

        # region page_home
        self.page_home = QWidget(self.stackedWidget)
        self.stackedWidget.addWidget(self.page_home)

        self.labelHome = QLabel(self.page_home)
        self.labelHome.setGeometry(QRect(180, 60, 440, 200))
        self.labelHome.setStyleSheet('font: 10pt "Consolas";')
        self.labelHome.setText(HOME)

        self.verticalLayoutWidget_home = QWidget(self.page_home)
        self.verticalLayoutWidget_home.setGeometry(QRect(320, 300, 160, 200))
        self.verticalLayout_home = QVBoxLayout(self.verticalLayoutWidget_home)
        self.verticalLayout_home.setContentsMargins(0, 0, 0, 0)

        self.pushButtonToInsert = QPushButton(self.verticalLayoutWidget_home)
        self.pushButtonToInsert.setObjectName("pushButtonToInsert")
        self.verticalLayout_home.addWidget(self.pushButtonToInsert)
        self.pushButtonToInsert.setText("Record insertion")

        self.pushButtonToUpdate = QPushButton(self.verticalLayoutWidget_home)
        self.pushButtonToUpdate.setObjectName("pushButtonToUpdate")
        self.verticalLayout_home.addWidget(self.pushButtonToUpdate)
        self.pushButtonToUpdate.setText("Record update")

        self.pushButtonToDelete = QPushButton(self.verticalLayoutWidget_home)
        self.pushButtonToDelete.setObjectName("pushButtonToDelete")
        self.verticalLayout_home.addWidget(self.pushButtonToDelete)
        self.pushButtonToDelete.setText("Record deletion")

        self.pushButtonToSearch = QPushButton(self.verticalLayoutWidget_home)
        self.pushButtonToSearch.setObjectName("pushButtonToSearch")
        self.verticalLayout_home.addWidget(self.pushButtonToSearch)
        self.pushButtonToSearch.setText("Record search")

        self.pushButtonToDistance = QPushButton(self.verticalLayoutWidget_home)
        self.pushButtonToDistance.setObjectName("pushButtonToDistance")
        self.verticalLayout_home.addWidget(self.pushButtonToDistance)
        self.pushButtonToDistance.setText("Genetic distance")
        # endregion

        # region page_insert
        self.page_insert = QWidget(self.stackedWidget)
        self.stackedWidget.addWidget(self.page_insert)

        self.comboBox_insert = QComboBox(self.page_insert)
        self.comboBox_insert.setObjectName("comboBox_insert")
        self.comboBox_insert.setGeometry(QRect(300, 60, 200, 24))
        self.comboBox_insert.addItem("Choose inserting data...")
        self.comboBox_insert.addItem("Insert SNP")
        self.comboBox_insert.addItem("Insert gene")
        self.comboBox_insert.addItem("Insert disease")

        self.line_insert = QFrame(self.page_insert)
        self.line_insert.setGeometry(QRect(200, 100, 400, 50))
        self.line_insert.setFrameShape(QFrame.HLine)
        self.line_insert.setFrameShadow(QFrame.Plain)

        self.stackedWidget_insert = QStackedWidget(self.page_insert)
        self.stackedWidget_insert.setGeometry(QRect(0, 150, 800, 450))

        # region page_insert_selectType
        self.page_insert_selectType = QWidget(self.stackedWidget_insert)
        self.stackedWidget_insert.addWidget(self.page_insert_selectType)

        self.label_insert = QLabel(self.page_insert_selectType)
        self.label_insert.setGeometry(QRect(180, 0, 440, 200))
        self.label_insert.setStyleSheet('font: 10pt "Consolas";')
        self.label_insert.setText(INSERT)
        # endregion

        # region page_insert_snp
        self.page_insert_snp = QWidget(self.stackedWidget_insert)
        self.stackedWidget_insert.addWidget(self.page_insert_snp)

        self.widget_insert_snp = InsertSnpWidget(self.page_insert_snp)
        self.widget_insert_snp.setGeometry(QRect(250, 25, 300, 300))

        self.pushButton_insert_snp = QPushButton(self.page_insert_snp)
        self.pushButton_insert_snp.setObjectName("pushButton_insert_snp")
        self.pushButton_insert_snp.setGeometry(350, 300, 100, 24)
        self.pushButton_insert_snp.setText("Submit")
        # endregion

        # region page_insert_gene
        self.page_insert_gene = QWidget(self.stackedWidget_insert)
        self.stackedWidget_insert.addWidget(self.page_insert_gene)

        self.widget_insert_gene = InsertGeneWidget(self.page_insert_gene)
        self.widget_insert_gene.setGeometry(QRect(250, 25, 300, 300))

        self.pushButton_insert_gene = QPushButton(self.page_insert_gene)
        self.pushButton_insert_gene.setObjectName("pushButton_insert_gene")
        self.pushButton_insert_gene.setGeometry(350, 300, 100, 24)
        self.pushButton_insert_gene.setText("Submit")
        # endregion

        # region page_insert_disease
        self.page_insert_disease = QWidget(self.stackedWidget_insert)
        self.stackedWidget_insert.addWidget(self.page_insert_disease)

        self.widget_insert_disease = InsertDiseaseWidget(self.page_insert_disease)
        self.widget_insert_disease.setGeometry(QRect(250, 25, 300, 300))

        self.pushButton_insert_disease = QPushButton(self.page_insert_disease)
        self.pushButton_insert_disease.setObjectName("pushButton_insert_disease")
        self.pushButton_insert_disease.setGeometry(350, 300, 100, 24)
        self.pushButton_insert_disease.setText("Submit")
        # endregion

        # endregion

        # region page_update
        self.page_update = QWidget(self.stackedWidget)
        self.stackedWidget.addWidget(self.page_update)

        self.comboBox_update = QComboBox(self.page_update)
        self.comboBox_update.setObjectName("comboBox_update")
        self.comboBox_update.setGeometry(QRect(300, 60, 200, 24))
        self.comboBox_update.addItem("Choose updating data...")
        self.comboBox_update.addItem("Update SNP")
        self.comboBox_update.addItem("Update gene")
        self.comboBox_update.addItem("Update disease")

        self.line_update = QFrame(self.page_update)
        self.line_update.setGeometry(QRect(200, 100, 400, 50))
        self.line_update.setFrameShape(QFrame.HLine)
        self.line_update.setFrameShadow(QFrame.Plain)

        self.stackedWidget_update = QStackedWidget(self.page_update)
        self.stackedWidget_update.setGeometry(QRect(0, 150, 800, 450))

        # region page_update_selectType
        self.page_update_selectType = QWidget(self.stackedWidget_update)
        self.stackedWidget_update.addWidget(self.page_update_selectType)

        self.label_update = QLabel(self.page_update_selectType)
        self.label_update.setGeometry(QRect(170, 0, 480, 200))
        self.label_update.setStyleSheet('font: 10pt "Consolas";')
        self.label_update.setText(UPDATE)
        # endregion

        # region page_update_snp
        self.page_update_snp = QWidget(self.stackedWidget_update)
        self.stackedWidget_update.addWidget(self.page_update_snp)

        self.label_update_snp = QLabel(self.page_update_snp)
        self.label_update_snp.setGeometry(QRect(200, 0, 400, 50))
        self.label_update_snp.setAlignment(Qt.AlignCenter)
        self.label_update_snp.setText("Type SNP ID below.")

        self.lineEdit_update_snp = QLineEdit(self.page_update_snp)
        self.lineEdit_update_snp.setGeometry(QRect(200, 50, 400, 25))
        self.lineEdit_update_snp.setValidator(QIntValidator(1, INT_MAX))
        self.lineEdit_update_snp.setAlignment(Qt.AlignCenter)
        self.lineEdit_update_snp.setPlaceholderText("positive integer")

        self.pushButton_update_snp = QPushButton(self.page_update_snp)
        self.pushButton_update_snp.setObjectName("pushButton_update_snp")
        self.pushButton_update_snp.setGeometry(350, 300, 100, 24)
        self.pushButton_update_snp.setText("Check")
        # endregion

        # region page_update_gene
        self.page_update_gene = QWidget(self.stackedWidget_update)
        self.stackedWidget_update.addWidget(self.page_update_gene)

        self.label_update_gene = QLabel(self.page_update_gene)
        self.label_update_gene.setGeometry(QRect(200, 0, 400, 50))
        self.label_update_gene.setAlignment(Qt.AlignCenter)
        self.label_update_gene.setText("Type Gene ID below.")

        self.lineEdit_update_gene = QLineEdit(self.page_update_gene)
        self.lineEdit_update_gene.setGeometry(QRect(200, 50, 400, 25))
        self.lineEdit_update_gene.setValidator(QIntValidator(1, INT_MAX))
        self.lineEdit_update_gene.setAlignment(Qt.AlignCenter)
        self.lineEdit_update_gene.setPlaceholderText("positive integer")

        self.pushButton_update_gene = QPushButton(self.page_update_gene)
        self.pushButton_update_gene.setObjectName("pushButton_update_gene")
        self.pushButton_update_gene.setGeometry(350, 300, 100, 24)
        self.pushButton_update_gene.setText("Check")
        # endregion

        # region page_update_disease
        self.page_update_disease = QWidget(self.stackedWidget_update)
        self.stackedWidget_update.addWidget(self.page_update_disease)

        self.label_update_disease = QLabel(self.page_update_disease)
        self.label_update_disease.setGeometry(QRect(200, 0, 400, 50))
        self.label_update_disease.setAlignment(Qt.AlignCenter)
        self.label_update_disease.setText("Type Disease ID below.")

        self.lineEdit_update_disease = QLineEdit(self.page_update_disease)
        self.lineEdit_update_disease.setGeometry(QRect(200, 50, 400, 25))
        self.lineEdit_update_disease.setValidator(QIntValidator(1, INT_MAX))
        self.lineEdit_update_disease.setAlignment(Qt.AlignCenter)
        self.lineEdit_update_disease.setPlaceholderText("positive integer")

        self.pushButton_update_disease = QPushButton(self.page_update_disease)
        self.pushButton_update_disease.setObjectName("pushButton_update_disease")
        self.pushButton_update_disease.setGeometry(350, 300, 100, 24)
        self.pushButton_update_disease.setText("Check")
        # endregion

        # region page_update_update
        self.page_update_update = QWidget(self.stackedWidget_update)
        self.stackedWidget_update.addWidget(self.page_update_update)
        self.widget_update_update = None  # used later
        self.pushButton_update_update = None  # used later
        # endregion

        # endregion

        # region page_delete
        self.page_delete = QWidget(self.stackedWidget)
        self.stackedWidget.addWidget(self.page_delete)

        self.comboBox_delete = QComboBox(self.page_delete)
        self.comboBox_delete.setObjectName("comboBox_delete")
        self.comboBox_delete.setGeometry(QRect(300, 60, 200, 24))
        self.comboBox_delete.addItem("Choose deleting data...")
        self.comboBox_delete.addItem("Delete SNP")
        self.comboBox_delete.addItem("Delete gene")
        self.comboBox_delete.addItem("Delete disease")

        self.line_delete = QFrame(self.page_delete)
        self.line_delete.setGeometry(QRect(200, 100, 400, 50))
        self.line_delete.setFrameShape(QFrame.HLine)
        self.line_delete.setFrameShadow(QFrame.Plain)

        self.stackedWidget_delete = QStackedWidget(self.page_delete)
        self.stackedWidget_delete.setGeometry(QRect(0, 150, 800, 450))

        # region page_delete_selectType
        self.page_delete_selectType = QWidget(self.stackedWidget_delete)
        self.stackedWidget_delete.addWidget(self.page_delete_selectType)

        self.label_delete = QLabel(self.page_delete_selectType)
        self.label_delete.setGeometry(QRect(190, 0, 440, 200))
        self.label_delete.setStyleSheet('font: 10pt "Consolas";')
        self.label_delete.setText(DELETE)
        # endregion

        # region page_delete_snp
        self.page_delete_snp = QWidget(self.stackedWidget_delete)
        self.stackedWidget_delete.addWidget(self.page_delete_snp)

        self.label_delete_snp = QLabel(self.page_delete_snp)
        self.label_delete_snp.setGeometry(QRect(200, 0, 400, 50))
        self.label_delete_snp.setAlignment(Qt.AlignCenter)
        self.label_delete_snp.setText("Type SNP ID below.")

        self.lineEdit_delete_snp = QLineEdit(self.page_delete_snp)
        self.lineEdit_delete_snp.setGeometry(QRect(200, 50, 400, 25))
        self.lineEdit_delete_snp.setValidator(QIntValidator(1, INT_MAX))
        self.lineEdit_delete_snp.setAlignment(Qt.AlignCenter)
        self.lineEdit_delete_snp.setPlaceholderText("positive integer")

        self.pushButton_delete_snp = QPushButton(self.page_delete_snp)
        self.pushButton_delete_snp.setObjectName("pushButton_delete_snp")
        self.pushButton_delete_snp.setGeometry(350, 300, 100, 24)
        self.pushButton_delete_snp.setText("Submit")
        # endregion

        # region page_delete_gene
        self.page_delete_gene = QWidget(self.stackedWidget_delete)
        self.stackedWidget_delete.addWidget(self.page_delete_gene)

        self.label_delete_gene = QLabel(self.page_delete_gene)
        self.label_delete_gene.setGeometry(QRect(200, 0, 400, 50))
        self.label_delete_gene.setAlignment(Qt.AlignCenter)
        self.label_delete_gene.setText("Type Gene ID below.")

        self.lineEdit_delete_gene = QLineEdit(self.page_delete_gene)
        self.lineEdit_delete_gene.setGeometry(QRect(200, 50, 400, 25))
        self.lineEdit_delete_gene.setValidator(QIntValidator(1, INT_MAX))
        self.lineEdit_delete_gene.setAlignment(Qt.AlignCenter)
        self.lineEdit_delete_gene.setPlaceholderText("positive integer")

        self.pushButton_delete_gene = QPushButton(self.page_delete_gene)
        self.pushButton_delete_gene.setObjectName("pushButton_delete_gene")
        self.pushButton_delete_gene.setGeometry(350, 300, 100, 24)
        self.pushButton_delete_gene.setText("Submit")
        # endregion

        # region page_delete_disease
        self.page_delete_disease = QWidget(self.stackedWidget_delete)
        self.stackedWidget_delete.addWidget(self.page_delete_disease)

        self.label_delete_disease = QLabel(self.page_delete_disease)
        self.label_delete_disease.setGeometry(QRect(200, 0, 400, 50))
        self.label_delete_disease.setAlignment(Qt.AlignCenter)
        self.label_delete_disease.setText("Type Disease ID below.")

        self.lineEdit_delete_disease = QLineEdit(self.page_delete_disease)
        self.lineEdit_delete_disease.setGeometry(QRect(200, 50, 400, 25))
        self.lineEdit_delete_disease.setValidator(QIntValidator(1, INT_MAX))
        self.lineEdit_delete_disease.setAlignment(Qt.AlignCenter)
        self.lineEdit_delete_disease.setPlaceholderText("positive integer")

        self.pushButton_delete_disease = QPushButton(self.page_delete_disease)
        self.pushButton_delete_disease.setObjectName("pushButton_delete_disease")
        self.pushButton_delete_disease.setGeometry(350, 300, 100, 24)
        self.pushButton_delete_disease.setText("Submit")
        # endregion

        # endregion

        # region page_search
        self.page_search = QWidget(self.stackedWidget)
        self.stackedWidget.addWidget(self.page_search)

        self.comboBox_search = QComboBox(self.page_search)
        self.comboBox_search.setObjectName("comboBox_search")
        self.comboBox_search.setGeometry(QRect(250, 60, 300, 24))
        self.comboBox_search.addItem("Choose searching type...")
        self.comboBox_search.addItem("Search in SNP table")
        self.comboBox_search.addItem("Search in gene table")
        self.comboBox_search.addItem("Search in disease table")
        self.comboBox_search.addItem("Search diseases by SNP ID")
        self.comboBox_search.addItem("Search SNP ID by disease name")

        self.line_search = QFrame(self.page_search)
        self.line_search.setGeometry(QRect(200, 100, 400, 50))
        self.line_search.setFrameShape(QFrame.HLine)
        self.line_search.setFrameShadow(QFrame.Plain)

        self.stackedWidget_search = QStackedWidget(self.page_search)
        self.stackedWidget_search.setGeometry(QRect(0, 150, 800, 450))

        # region page_search_selectType
        self.page_search_selectType = QWidget(self.stackedWidget_search)
        self.stackedWidget_search.addWidget(self.page_search_selectType)

        self.label_search = QLabel(self.page_search_selectType)
        self.label_search.setGeometry(QRect(150, 0, 500, 200))
        self.label_search.setStyleSheet('font: 10pt "Consolas";')
        self.label_search.setText(SEARCH)
        # endregion

        # region page_search_snp
        self.page_search_snp = QWidget(self.stackedWidget_search)
        self.stackedWidget_search.addWidget(self.page_search_snp)

        self.widget_search_snp = InsertSnpWidget(self.page_search_snp)
        self.widget_search_snp.setGeometry(QRect(250, 25, 300, 300))
        # add exceptional cases
        self.widget_search_snp.comboBox_insert_snp_anc.insertItem(0, "")
        self.widget_search_snp.comboBox_insert_snp_anc.setCurrentIndex(0)
        self.widget_search_snp.comboBox_insert_snp_minor.insertItem(0, "")
        self.widget_search_snp.comboBox_insert_snp_minor.setCurrentIndex(0)

        self.label_search_snp = QLabel(self.page_search_snp)
        self.label_search_snp.setGeometry(QRect(200, 250, 400, 50))
        self.label_search_snp.setAlignment(Qt.AlignCenter)
        self.label_search_snp.setText("Leave empty to allow any values for the field.")

        self.pushButton_search_snp = QPushButton(self.page_search_snp)
        self.pushButton_search_snp.setObjectName("pushButton_search_snp")
        self.pushButton_search_snp.setGeometry(350, 300, 100, 24)
        self.pushButton_search_snp.setText("Submit")
        # endregion

        # region page_search_gene
        self.page_search_gene = QWidget(self.stackedWidget_search)
        self.stackedWidget_search.addWidget(self.page_search_gene)

        self.widget_search_gene = InsertGeneWidget(self.page_search_gene)
        self.widget_search_gene.setGeometry(QRect(250, 25, 300, 300))
        # add exceptional cases
        self.widget_search_gene.form.removeRow(8)

        self.label_search_gene = QLabel(self.page_search_gene)
        self.label_search_gene.setGeometry(QRect(200, 250, 400, 50))
        self.label_search_gene.setAlignment(Qt.AlignCenter)
        self.label_search_gene.setText("Leave empty to allow any values for the field.")

        self.pushButton_search_gene = QPushButton(self.page_search_gene)
        self.pushButton_search_gene.setObjectName("pushButton_search_gene")
        self.pushButton_search_gene.setGeometry(350, 300, 100, 24)
        self.pushButton_search_gene.setText("Submit")
        # endregion

        # region page_search_disease
        self.page_search_disease = QWidget(self.stackedWidget_search)
        self.stackedWidget_search.addWidget(self.page_search_disease)

        self.widget_search_disease = InsertDiseaseWidget(self.page_search_disease)
        self.widget_search_disease.setGeometry(QRect(250, 25, 300, 300))

        self.label_search_disease = QLabel(self.page_search_disease)
        self.label_search_disease.setGeometry(QRect(200, 250, 400, 50))
        self.label_search_disease.setAlignment(Qt.AlignCenter)
        self.label_search_disease.setText(
            "Leave empty to allow any values for the field."
        )

        self.pushButton_search_disease = QPushButton(self.page_search_disease)
        self.pushButton_search_disease.setObjectName("pushButton_search_disease")
        self.pushButton_search_disease.setGeometry(350, 300, 100, 24)
        self.pushButton_search_disease.setText("Submit")
        # endregion

        # region page_search_diseaseBySnpID
        self.page_search_diseaseBySnpID = QWidget(self.stackedWidget_search)
        self.stackedWidget_search.addWidget(self.page_search_diseaseBySnpID)

        self.label_search_diseaseBySnpID = QLabel(self.page_search_diseaseBySnpID)
        self.label_search_diseaseBySnpID.setGeometry(QRect(200, 0, 400, 50))
        self.label_search_diseaseBySnpID.setAlignment(Qt.AlignCenter)
        self.label_search_diseaseBySnpID.setText("Type SNP ID below.")

        self.lineEdit_search_diseaseBySnpID = QLineEdit(self.page_search_diseaseBySnpID)
        self.lineEdit_search_diseaseBySnpID.setGeometry(QRect(200, 50, 400, 25))
        self.lineEdit_search_diseaseBySnpID.setValidator(QIntValidator(1, INT_MAX))
        self.lineEdit_search_diseaseBySnpID.setAlignment(Qt.AlignCenter)
        self.lineEdit_search_diseaseBySnpID.setPlaceholderText("positive integer")

        self.pushButton_search_diseaseBySnpID = QPushButton(
            self.page_search_diseaseBySnpID
        )
        self.pushButton_search_diseaseBySnpID.setObjectName(
            "pushButton_search_diseaseBySnpID"
        )
        self.pushButton_search_diseaseBySnpID.setGeometry(350, 300, 100, 24)
        self.pushButton_search_diseaseBySnpID.setText("Submit")
        # endregion

        # region page_search_snpByDiseaseName
        self.page_search_snpByDiseaseName = QWidget(self.stackedWidget_search)
        self.stackedWidget_search.addWidget(self.page_search_snpByDiseaseName)

        self.label_search_snpByDiseaseName = QLabel(self.page_search_snpByDiseaseName)
        self.label_search_snpByDiseaseName.setGeometry(QRect(200, 0, 400, 50))
        self.label_search_snpByDiseaseName.setAlignment(Qt.AlignCenter)
        self.label_search_snpByDiseaseName.setText("Type disease name below.")

        self.lineEdit_search_snpByDiseaseName = QLineEdit(
            self.page_search_snpByDiseaseName
        )
        self.lineEdit_search_snpByDiseaseName.setGeometry(QRect(200, 50, 400, 25))
        self.lineEdit_search_snpByDiseaseName.setAlignment(Qt.AlignCenter)
        self.lineEdit_search_snpByDiseaseName.setPlaceholderText(
            "disease name (string)"
        )

        self.pushButton_search_snpByDiseaseName = QPushButton(
            self.page_search_snpByDiseaseName
        )
        self.pushButton_search_snpByDiseaseName.setObjectName(
            "pushButton_search_snpByDiseaseName"
        )
        self.pushButton_search_snpByDiseaseName.setGeometry(350, 300, 100, 24)
        self.pushButton_search_snpByDiseaseName.setText("Submit")
        # endregion

        # region page_search_result
        self.page_search_result = QWidget(self.stackedWidget_search)
        self.stackedWidget_search.addWidget(self.page_search_result)

        self.table_search_result = QTableWidget(self.page_search_result)
        self.table_search_result.setGeometry(QRect(100, 25, 600, 325))
        self.table_search_result.setEditTriggers(QAbstractItemView.NoEditTriggers)

        # overlapped
        self.widget_search_result_figure = QWidget(self.page_search_result)
        self.widget_search_result_figure.setGeometry(QRect(198, 275, 404, 96))

        self.figure_search_result = None  # used later
        # endregion

        # endregion

        # region page_distance
        self.page_distance = QWidget(self.stackedWidget)
        self.stackedWidget.addWidget(self.page_distance)

        self.radioButton_gds = QRadioButton(
            "Genetic distances among counties", self.page_distance
        )
        self.radioButton_gds.move(200, 60)
        self.radioButton_gds.setChecked(True)
        self.radioButton_gds.clicked.connect(self.on_radioButton_distance_clicked)

        self.radioButton_ind = QRadioButton(
            "Individual investigation", self.page_distance
        )
        self.radioButton_ind.move(450, 60)
        self.radioButton_ind.clicked.connect(self.on_radioButton_distance_clicked)

        self.line_distance = QFrame(self.page_distance)
        self.line_distance.setGeometry(QRect(200, 100, 400, 50))
        self.line_distance.setFrameShape(QFrame.HLine)
        self.line_distance.setFrameShadow(QFrame.Plain)

        self.stackedWidget_distance = QStackedWidget(self.page_distance)
        self.stackedWidget_distance.setGeometry(QRect(0, 150, 800, 450))

        # region page_distance_gds
        self.page_distance_gds = QWidget(self.stackedWidget_distance)
        self.stackedWidget_distance.addWidget(self.page_distance_gds)

        self.comboBox_distance_gds = QComboBox(self.page_distance_gds)
        self.comboBox_distance_gds.setObjectName("comboBox_distance_gds")
        self.comboBox_distance_gds.setGeometry(QRect(300, 0, 200, 24))
        self.comboBox_distance_gds.addItem("Choose a county...")
        county_list = [
            "Bartow",
            "Sussex",
            "KingandQueen",
            "NewKent",
            "Warren",
            "Northampton",
            "Columbus",
            "Williamsburg",
            "Onslow",
            "Georgetown",
            "Beaufort",
            "Bertie",
            "Craven",
            "Chatham",
            "Carteret",
            "Tyrrell",
            "Brunswick",
            "Martin",
            "Hertford",
            "Chowan",
            "Greene",
            "Tuscaloosa",
            "Hale",
            "Bibb",
            "Pickens",
            "Sumter",
            "Anson",
            "Jasper",
            "Tallapoosa",
            "PrinceGeorge",
            "Dinwiddie",
            "Mccurtain",
            "Cleveland",
            "Jones",
            "Franklin",
            "Marshall",
            "Winston",
            "Lamar",
            "Montgomery",
            "Richmond",
            "Hancock",
            "Fayette",
            "Dorchester",
            "Berkeley",
            "Saluda",
            "Moore",
            "Newberry",
            "Wake",
            "Lee",
            "Marengo",
            "Choctaw",
            "Clarke",
            "Greenwood",
            "Marion",
            "Chesterfield",
            "Levy",
            "Escambia",
            "Conecuh",
            "Lunenburg",
            "Mecklenburg",
            "Cumberland",
            "PrinceEdward",
            "Halifax",
            "Rowan",
            "Durham",
            "Gates",
            "Pasquotank",
            "Bamberg",
            "Wilcox",
            "Dallas",
            "Heard",
            "Calhoun",
            "Union",
            "Aiken",
            "York",
            "Edgefield",
            "Chattooga",
            "Spotsylvania",
            "Fluvanna",
            "Houston",
            "Polk",
            "Tyler",
            "Morehouse",
            "Pike",
            "Horry",
        ]
        for cou in county_list:
            self.comboBox_distance_gds.addItem(cou)

        self.stackedWidget_distance_gds = QStackedWidget(self.page_distance_gds)
        self.stackedWidget_distance_gds.setGeometry(QRect(0, 150, 800, 350))
        # endregion

        # region page_distance_gds_plot
        self.page_distance_gds_plot = QWidget(self.stackedWidget_distance_gds)
        self.stackedWidget_distance_gds.addWidget(self.page_distance_gds_plot)

        self.pushButton_distance_gds_plot = QPushButton(self.page_distance_gds_plot)
        self.pushButton_distance_gds_plot.setObjectName("pushButton_distance_gds_plot")
        self.pushButton_distance_gds_plot.setGeometry(300, 0, 200, 24)
        self.pushButton_distance_gds_plot.setText("Show Figure")
        # endregion

        # region page_distance_ind
        self.page_distance_ind = QWidget(self.stackedWidget_distance)
        self.stackedWidget_distance.addWidget(self.page_distance_ind)

        self.lineEdit_distance_ind = QLineEdit(self.page_distance_ind)
        self.lineEdit_distance_ind.setGeometry(QRect(200, 100, 400, 25))
        self.lineEdit_distance_ind.setAlignment(Qt.AlignCenter)
        self.lineEdit_distance_ind.setPlaceholderText("SNP genotype input file name")

        self.pushButton_distance_ind = QPushButton(self.page_distance_ind)
        self.pushButton_distance_ind.setObjectName("pushButton_distance_ind")
        self.pushButton_distance_ind.setGeometry(350, 300, 100, 24)
        self.pushButton_distance_ind.setText("Submit")
        # endregion

        # region page_distance_ind_result
        self.page_distance_ind_result = QWidget(self.stackedWidget_distance)
        self.stackedWidget_distance.addWidget(self.page_distance_ind_result)

        self.pushButton_distance_ind_plot = QPushButton(self.page_distance_ind_result)
        self.pushButton_distance_ind_plot.setObjectName("pushButton_distance_ind_plot")
        self.pushButton_distance_ind_plot.setGeometry(300, 0, 200, 24)
        self.pushButton_distance_ind_plot.setText("Show Figure")
        # endregion

        # endregion

        QMetaObject.connectSlotsByName(self)

    # region connected functions

    # region action
    @Slot()
    def on_actionHome_triggered(self):
        self.stackedWidget.setCurrentWidget(self.page_home)

    @Slot()
    def on_actionOpen_triggered(self):
        MainWindow().show()

    @Slot()
    def on_actionClose_triggered(self):
        self.close()

    @Slot()
    def on_actionShowLog_triggered(self):
        global log_window
        if log_window is None:
            LogWindow().show()
        else:
            assert isinstance(log_window, LogWindow)
            log_window.setWindowState(Qt.WindowActive)

    @Slot()
    def on_actionExportLog_triggered(self):
        exportLog()

    @Slot()
    def on_actionSQL_triggered(self):
        sql_window = SQLwindow()
        sql_window.exec()

    # endregion

    # region page_home
    @Slot()
    def on_pushButtonToInsert_clicked(self):
        self.stackedWidget.setCurrentWidget(self.page_insert)
        self.comboBox_insert.setCurrentIndex(0)
        self.stackedWidget_insert.setCurrentIndex(0)

    @Slot()
    def on_pushButtonToUpdate_clicked(self):
        self.stackedWidget.setCurrentWidget(self.page_update)
        self.comboBox_update.setCurrentIndex(0)
        self.stackedWidget_update.setCurrentIndex(0)

    @Slot()
    def on_pushButtonToDelete_clicked(self):
        self.stackedWidget.setCurrentWidget(self.page_delete)
        self.comboBox_delete.setCurrentIndex(0)
        self.stackedWidget_delete.setCurrentIndex(0)

    @Slot()
    def on_pushButtonToSearch_clicked(self):
        self.stackedWidget.setCurrentWidget(self.page_search)
        self.comboBox_search.setCurrentIndex(0)
        self.stackedWidget_search.setCurrentIndex(0)

    @Slot()
    def on_pushButtonToDistance_clicked(self):
        self.stackedWidget.setCurrentWidget(self.page_distance)

    # endregion

    # region page_insert
    @Slot(int)
    def on_comboBox_insert_currentIndexChanged(self, index):
        self.stackedWidget_insert.setCurrentIndex(index)

    @Slot()
    def on_pushButton_insert_snp_clicked(self):
        self.statusbar.showMessage("Validating fields...")
        if self.formCheckValid(self.widget_insert_snp.form):
            if (
                self.widget_insert_snp.comboBox_insert_snp_anc.currentText()
                == self.widget_insert_snp.comboBox_insert_snp_minor.currentText()
            ):
                self.statusbar.showMessage("Invalid field: alleles", TIMEOUT_MS)
                return
            self.statusbar.showMessage("Parsing...")
            value_tuple = formToTuple(self.widget_insert_snp.form)
            query = "INSERT INTO snp VALUES" + str(value_tuple) + ";"
            self.statusbar.showMessage("Sending to server...")
            self.statusbar.showMessage(getResponse(query), TIMEOUT_MS)

    @Slot()
    def on_pushButton_insert_gene_clicked(self):
        self.statusbar.showMessage("Validating fields...")
        if self.formCheckValid(self.widget_insert_gene.form):
            self.statusbar.showMessage("Parsing...")
            value_tuple = formToTuple(self.widget_insert_gene.form)
            query = "INSERT INTO gene VALUES" + str(value_tuple) + ";"
            self.statusbar.showMessage("Sending to server...")
            self.statusbar.showMessage(getResponse(query), TIMEOUT_MS)

    @Slot()
    def on_pushButton_insert_disease_clicked(self):
        self.statusbar.showMessage("Validating fields...")
        if self.formCheckValid(self.widget_insert_disease.form):
            self.statusbar.showMessage("Parsing...")
            value_tuple = formToTuple(self.widget_insert_disease.form)
            query = "INSERT INTO disease VALUES" + str(value_tuple) + ";"
            self.statusbar.showMessage("Sending to server...")
            self.statusbar.showMessage(getResponse(query), TIMEOUT_MS)

    # endregion

    # region page_update
    @Slot(int)
    def on_comboBox_update_currentIndexChanged(self, index):
        self.clearPageUpdateUpdate()
        self.stackedWidget_update.setCurrentIndex(index)

    @Slot()
    def on_pushButton_update_snp_clicked(self):
        self.statusbar.showMessage("Validating fields...")
        if fieldCheckValid(self.lineEdit_update_snp):
            self.statusbar.showMessage("Parsing...")
            query = (
                "SELECT * FROM snp WHERE snp_id = "
                + self.lineEdit_update_snp.text()
                + ";"
            )
            self.statusbar.showMessage("Sending to server...")
            self.dataBuffer = getResponse(query)
            self.statusbar.showMessage("Showing search result...", TIMEOUT_MS)
            assert isinstance(self.dataBuffer, DataFrame)
            assert 0 <= len(self.dataBuffer) <= 1
            if len(self.dataBuffer) == 0:
                self.statusbar.showMessage("No such record.", TIMEOUT_MS)
            else:
                self.widget_update_update = InsertSnpWidget(self.page_update_update)
                self.widget_update_update.setGeometry(QRect(250, 25, 300, 300))

                self.pushButton_update_update = QPushButton(self.page_update_update)
                self.pushButton_update_update.setObjectName("pushButton_update_update")
                self.pushButton_update_update.setGeometry(350, 300, 100, 24)
                self.pushButton_update_update.setText("Update")
                self.pushButton_update_update.clicked.connect(
                    self.on_pushButton_update_update_snp_clicked
                )

                self.widget_update_update.lineEdit_insert_snp_id.setEnabled(False)
                searchResultToForm(self.widget_update_update.form, self.dataBuffer)
                self.stackedWidget_update.setCurrentWidget(self.page_update_update)
        else:
            self.statusbar.showMessage("Invalid field: SNP ID", TIMEOUT_MS)

    def on_pushButton_update_update_snp_clicked(self):
        assert isinstance(self.widget_update_update, InsertSnpWidget)
        self.statusbar.showMessage("Validating fields...")
        if self.formCheckValid(self.widget_update_update.form):
            if (
                self.widget_update_update.comboBox_insert_snp_anc.currentText()
                == self.widget_update_update.comboBox_insert_snp_minor.currentText()
            ):
                self.statusbar.showMessage("Invalid field: alleles", TIMEOUT_MS)
                return
            self.statusbar.showMessage("Parsing...")
            value_tuple = formToTuple(self.widget_update_update.form)
            wrapped_tuple = tuple(map(quoteWrappedStr, value_tuple))
            query = "UPDATE snp SET "
            for i in range(len(wrapped_tuple)):
                if i == PKIDX_SNP:
                    continue
                query = "".join((query, FIELD_SNP[i], "=", wrapped_tuple[i], ", "))
            query = "".join(
                (
                    query[:-2],
                    " WHERE ",
                    FIELD_SNP[PKIDX_SNP],
                    "=",
                    wrapped_tuple[PKIDX_SNP],
                    ";",
                )
            )
            self.statusbar.showMessage("Sending to server...")
            self.statusbar.showMessage(getResponse(query), TIMEOUT_MS)
            self.stackedWidget_update.setCurrentWidget(self.page_update_snp)
            self.clearPageUpdateUpdate()

    @Slot()
    def on_pushButton_update_gene_clicked(self):
        self.statusbar.showMessage("Validating fields...")
        if fieldCheckValid(self.lineEdit_update_gene):
            self.statusbar.showMessage("Parsing...")
            query = (
                "SELECT * FROM gene WHERE gene_id = "
                + self.lineEdit_update_gene.text()
                + ";"
            )
            self.statusbar.showMessage("Sending to server...")
            self.dataBuffer = getResponse(query)
            self.statusbar.showMessage("Showing search result...", TIMEOUT_MS)
            assert isinstance(self.dataBuffer, DataFrame)
            assert 0 <= len(self.dataBuffer) <= 1
            if len(self.dataBuffer) == 0:
                self.statusbar.showMessage("No such record.", TIMEOUT_MS)
            else:
                self.widget_update_update = InsertGeneWidget(self.page_update_update)
                self.widget_update_update.setGeometry(QRect(250, 25, 300, 300))

                self.pushButton_update_update = QPushButton(self.page_update_update)
                self.pushButton_update_update.setObjectName("pushButton_update_update")
                self.pushButton_update_update.setGeometry(350, 300, 100, 24)
                self.pushButton_update_update.setText("Update")
                self.pushButton_update_update.clicked.connect(
                    self.on_pushButton_update_update_gene_clicked
                )

                self.widget_update_update.lineEdit_insert_gene_id.setEnabled(False)
                searchResultToForm(self.widget_update_update.form, self.dataBuffer)
                self.stackedWidget_update.setCurrentWidget(self.page_update_update)
        else:
            self.statusbar.showMessage("Invalid field: Gene ID", TIMEOUT_MS)

    def on_pushButton_update_update_gene_clicked(self):
        assert isinstance(self.widget_update_update, InsertGeneWidget)
        self.statusbar.showMessage("Validating fields...")
        if self.formCheckValid(self.widget_update_update.form):
            self.statusbar.showMessage("Parsing...")
            value_tuple = formToTuple(self.widget_update_update.form)
            wrapped_tuple = tuple(map(quoteWrappedStr, value_tuple))
            query = "UPDATE gene SET "
            for i in range(len(wrapped_tuple)):
                if i == PKIDX_GENE:
                    continue
                query = "".join((query, FIELD_GENE[i], "=", wrapped_tuple[i], ", "))
            query = "".join(
                (
                    query[:-2],
                    " WHERE ",
                    FIELD_GENE[PKIDX_GENE],
                    "=",
                    wrapped_tuple[PKIDX_GENE],
                    ";",
                )
            )
            self.statusbar.showMessage("Sending to server...")
            self.statusbar.showMessage(getResponse(query), TIMEOUT_MS)
            self.stackedWidget_update.setCurrentWidget(self.page_update_gene)
            self.clearPageUpdateUpdate()

    @Slot()
    def on_pushButton_update_disease_clicked(self):
        self.statusbar.showMessage("Validating fields...")
        if fieldCheckValid(self.lineEdit_update_disease):
            self.statusbar.showMessage("Parsing...")
            query = (
                "SELECT * FROM disease WHERE disease_id = "
                + self.lineEdit_update_disease.text()
                + ";"
            )
            self.statusbar.showMessage("Sending to server...")
            self.dataBuffer = getResponse(query)
            self.statusbar.showMessage("Showing search result...", TIMEOUT_MS)
            assert isinstance(self.dataBuffer, DataFrame)
            assert 0 <= len(self.dataBuffer) <= 1
            if len(self.dataBuffer) == 0:
                self.statusbar.showMessage("No such record.", TIMEOUT_MS)
            else:
                self.widget_update_update = InsertDiseaseWidget(self.page_update_update)
                self.widget_update_update.setGeometry(QRect(250, 25, 300, 300))

                self.pushButton_update_update = QPushButton(self.page_update_update)
                self.pushButton_update_update.setObjectName("pushButton_update_update")
                self.pushButton_update_update.setGeometry(350, 300, 100, 24)
                self.pushButton_update_update.setText("Update")
                self.pushButton_update_update.clicked.connect(
                    self.on_pushButton_update_update_disease_clicked
                )

                self.widget_update_update.lineEdit_insert_disease_id.setEnabled(False)
                searchResultToForm(self.widget_update_update.form, self.dataBuffer)
                self.stackedWidget_update.setCurrentWidget(self.page_update_update)
        else:
            self.statusbar.showMessage("Invalid field: Disease ID", TIMEOUT_MS)

    def on_pushButton_update_update_disease_clicked(self):
        assert isinstance(self.widget_update_update, InsertDiseaseWidget)
        self.statusbar.showMessage("Validating fields...")
        if self.formCheckValid(self.widget_update_update.form):
            self.statusbar.showMessage("Parsing...")
            value_tuple = formToTuple(self.widget_update_update.form)
            wrapped_tuple = tuple(map(quoteWrappedStr, value_tuple))
            query = "UPDATE disease SET "
            for i in range(len(wrapped_tuple)):
                if i == PKIDX_DISEASE:
                    continue
                query = "".join((query, FIELD_DISEASE[i], "=", wrapped_tuple[i], ", "))
            query = "".join(
                (
                    query[:-2],
                    " WHERE ",
                    FIELD_DISEASE[PKIDX_DISEASE],
                    "=",
                    wrapped_tuple[PKIDX_DISEASE],
                    ";",
                )
            )
            self.statusbar.showMessage("Sending to server...")
            self.statusbar.showMessage(getResponse(query), TIMEOUT_MS)
            self.stackedWidget_update.setCurrentWidget(self.page_update_disease)
            self.clearPageUpdateUpdate()

    # endregion

    # region page_delete
    @Slot(int)
    def on_comboBox_delete_currentIndexChanged(self, index):
        self.stackedWidget_delete.setCurrentIndex(index)

    @Slot()
    def on_pushButton_delete_snp_clicked(self):
        self.statusbar.showMessage("Validating fields...")
        if fieldCheckValid(self.lineEdit_delete_snp):
            self.statusbar.showMessage("Parsing...")
            query = (
                "DELETE FROM snp WHERE snp_id = "
                + self.lineEdit_delete_snp.text()
                + ";"
            )
            self.statusbar.showMessage("Sending to server...")
            self.statusbar.showMessage(getResponse(query), TIMEOUT_MS)
        else:
            self.statusbar.showMessage("Invalid field: SNP ID", TIMEOUT_MS)

    @Slot()
    def on_pushButton_delete_gene_clicked(self):
        self.statusbar.showMessage("Validating fields...")
        if fieldCheckValid(self.lineEdit_delete_gene):
            self.statusbar.showMessage("Parsing...")
            query = (
                "DELETE FROM gene WHERE gene_id = "
                + self.lineEdit_delete_gene.text()
                + ";"
            )
            self.statusbar.showMessage("Sending to server...")
            self.statusbar.showMessage(getResponse(query), TIMEOUT_MS)
        else:
            self.statusbar.showMessage("Invalid field: Gene ID", TIMEOUT_MS)

    @Slot()
    def on_pushButton_delete_disease_clicked(self):
        self.statusbar.showMessage("Validating fields...")
        if fieldCheckValid(self.lineEdit_delete_disease):
            self.statusbar.showMessage("Parsing...")
            query = (
                "DELETE FROM disease WHERE disease_id = "
                + self.lineEdit_delete_disease.text()
                + ";"
            )
            self.statusbar.showMessage("Sending to server...")
            self.statusbar.showMessage(getResponse(query), TIMEOUT_MS)
        else:
            self.statusbar.showMessage("Invalid field: Disease ID", TIMEOUT_MS)

    # endregion

    # region page_search
    @Slot(int)
    def on_comboBox_search_currentIndexChanged(self, index):
        self.clearSearchResultSnpFigure()
        self.stackedWidget_search.setCurrentIndex(index)

    @Slot()
    def on_pushButton_search_snp_clicked(self):
        self.statusbar.showMessage("Parsing...")
        query = "SELECT * FROM snp WHERE "
        fields = formToTuple(self.widget_search_snp.form, True)
        flag = False
        for field in fields:
            if field != "":
                flag = True
                break
        if flag:
            for i in range(len(fields)):
                field = fields[i]
                if field == "":
                    continue
                query = "".join(
                    (query, FIELD_SNP[i], "=", quoteWrappedStr(field), " AND ")
                )
            query = query[:-5] + ";"
        else:
            query = "SELECT * FROM snp"

        self.statusbar.showMessage("Sending to server...")
        self.dataBuffer = getResponse(query)
        self.statusbar.showMessage("Showing search result...", TIMEOUT_MS)
        assert isinstance(self.dataBuffer, DataFrame)
        if len(self.dataBuffer) == 0:
            self.statusbar.showMessage("No such record.", TIMEOUT_MS)
        else:
            searchResultToTable(self.table_search_result, self.dataBuffer, LABEL_SNP)
            self.table_search_result.resize(600, 325)
            self.table_search_result.resizeColumnsToContents()
            try:
                self.table_search_result.cellDoubleClicked.disconnect()
            except Exception:
                pass
            self.stackedWidget_search.setCurrentWidget(self.page_search_result)

    @Slot()
    def on_pushButton_search_gene_clicked(self):
        self.statusbar.showMessage("Parsing...")
        query = "SELECT * FROM gene WHERE "
        fields = formToTuple(self.widget_search_gene.form, True)
        flag = False
        for field in fields:
            if field != "":
                flag = True
                break
        if flag:
            for i in range(len(fields)):
                field = fields[i]
                if field == "":
                    continue
                query = "".join(
                    (query, FIELD_GENE[i], "=", quoteWrappedStr(field), " AND ")
                )
            query = query[:-5] + ";"
        else:
            query = "SELECT * FROM gene"

        self.statusbar.showMessage("Sending to server...")
        self.dataBuffer = getResponse(query)
        self.statusbar.showMessage("Showing search result...", TIMEOUT_MS)
        assert isinstance(self.dataBuffer, DataFrame)
        if len(self.dataBuffer) == 0:
            self.statusbar.showMessage("No such record.", TIMEOUT_MS)
        else:
            searchResultToTable(self.table_search_result, self.dataBuffer, LABEL_GENE)
            self.table_search_result.resize(600, 325)
            self.table_search_result.resizeColumnsToContents()
            try:
                self.table_search_result.cellDoubleClicked.disconnect()
            except Exception:
                pass
            self.stackedWidget_search.setCurrentWidget(self.page_search_result)

    @Slot()
    def on_pushButton_search_disease_clicked(self):
        self.statusbar.showMessage("Parsing...")
        fields = formToTuple(self.widget_search_disease.form, True)
        if fields[0] == "" and fields[1] == "":
            query = "SELECT * FROM disease;"
        elif fields[0] == "" and fields[1] != "":
            query = "SELECT * FROM disease WHERE name LIKE '%%%s%%';" % fields[1]
        elif fields[0] != "" and fields[1] == "":
            query = "SELECT * FROM disease WHERE disease_id=%d;" % fields[0]
        else:
            query = (
                "SELECT * FROM disease WHERE disease_id=%d AND name LIKE '%%%s%%';"
                % (
                    fields[0],
                    fields[1],
                )
            )

        self.statusbar.showMessage("Sending to server...")
        self.dataBuffer = getResponse(query)
        self.statusbar.showMessage("Showing search result...", TIMEOUT_MS)
        assert isinstance(self.dataBuffer, DataFrame)
        if len(self.dataBuffer) == 0:
            self.statusbar.showMessage("No such record.", TIMEOUT_MS)
        else:
            searchResultToTable(
                self.table_search_result, self.dataBuffer, LABEL_DISEASE
            )
            self.table_search_result.resize(600, 325)
            self.table_search_result.resizeColumnsToContents()
            try:
                self.table_search_result.cellDoubleClicked.disconnect()
            except Exception:
                pass
            self.stackedWidget_search.setCurrentWidget(self.page_search_result)

    @Slot()
    def on_pushButton_search_diseaseBySnpID_clicked(self):
        if not fieldCheckValid(self.lineEdit_search_diseaseBySnpID):
            self.statusbar.showMessage("Invalid field: SNP ID", TIMEOUT_MS)
            return
        self.statusbar.showMessage("Parsing...")
        query = (
            "SELECT snp_id, disease_id, name FROM snp_gene \
                LEFT OUTER JOIN gene_disease USING (gene_id) \
                    RIGHT OUTER JOIN disease USING (disease_id) \
                        RIGHT OUTER JOIN snp USING (snp_id) WHERE snp_id=%s;"
            % self.lineEdit_search_diseaseBySnpID.text()
        )

        self.statusbar.showMessage("Sending to server...")
        self.dataBuffer = getResponse(query)
        self.statusbar.showMessage("Showing search result...", TIMEOUT_MS)
        assert isinstance(self.dataBuffer, DataFrame)
        if len(self.dataBuffer) == 0:
            self.statusbar.showMessage("No such record.", TIMEOUT_MS)
        else:
            searchResultToTable(
                self.table_search_result,
                self.dataBuffer,
                ("SNP ID", "Disease ID", "Name"),
            )

            # showing figure
            self.table_search_result.resize(600, 200)
            self.table_search_result.resizeColumnsToContents()
            self.statusbar.showMessage(
                "Double-click cell to view SNP figure.", TIMEOUT_MS
            )
            try:
                self.table_search_result.cellDoubleClicked.disconnect()
            except Exception:
                pass
            self.table_search_result.cellDoubleClicked.connect(
                self.on_table_search_result_cellDoubleClicked
            )
            self.stackedWidget_search.setCurrentWidget(self.page_search_result)

    @Slot()
    def on_pushButton_search_snpByDiseaseName_clicked(self):
        if self.lineEdit_search_snpByDiseaseName.text() == "":
            self.statusbar.showMessage("Empty keyword is not allowed.", TIMEOUT_MS)
            return
        self.statusbar.showMessage("Parsing...")
        query = (
            "SELECT snp_id, disease_id, name FROM snp_gene \
                RIGHT OUTER JOIN gene_disease USING (gene_id) \
                    RIGHT OUTER JOIN disease USING (disease_id) WHERE name LIKE '%%%s%%';"
            % self.lineEdit_search_snpByDiseaseName.text()
        )

        self.statusbar.showMessage("Sending to server...")
        self.dataBuffer = getResponse(query)
        self.statusbar.showMessage("Showing search result...", TIMEOUT_MS)
        assert isinstance(self.dataBuffer, DataFrame)
        if len(self.dataBuffer) == 0:
            self.statusbar.showMessage("No such record.", TIMEOUT_MS)
        else:
            searchResultToTable(
                self.table_search_result,
                self.dataBuffer,
                ("SNP ID", "Disease ID", "Name"),
            )

            # showing figure
            self.table_search_result.resize(600, 200)
            self.table_search_result.resizeColumnsToContents()
            self.statusbar.showMessage(
                "Double-click cell to view SNP figure.", TIMEOUT_MS
            )
            try:
                self.table_search_result.cellDoubleClicked.disconnect()
            except Exception:
                pass
            self.table_search_result.cellDoubleClicked.connect(
                self.on_table_search_result_cellDoubleClicked
            )
            self.stackedWidget_search.setCurrentWidget(self.page_search_result)

    def on_table_search_result_cellDoubleClicked(self, row, column):
        self.statusbar.showMessage("Showing SNP figure...", TIMEOUT_MS)
        snp = self.table_search_result.item(row, 0).text()
        self.clearSearchResultSnpFigure()
        try:
            snp_id = int(snp)
            self.figure_search_result = display_snp(snp_id)
            self.figure_search_result.setParent(self.widget_search_result_figure)
            self.figure_search_result.show()
        except ValueError:
            self.statusbar.showMessage("No corresponding SNP.", TIMEOUT_MS)

    # endregion

    # region page_distance
    def on_radioButton_distance_clicked(self):
        assert self.radioButton_gds.isChecked or self.radioButton_gds.isChecked()
        if self.radioButton_gds.isChecked():
            self.stackedWidget_distance.setCurrentIndex(0)
        else:
            self.stackedWidget_distance.setCurrentIndex(1)

    # endregion

    # region page_distance_gds
    @Slot(int)
    def on_comboBox_distance_gds_currentIndexChanged(self, index):
        cou = self.comboBox_distance_gds.currentIndex()
        if cou == 0:
            data = GD.create_county_dist_data(
                GD.county, GD.county_fips, GD.county_dist[0]
            )
        else:
            cou = cou - 1
            data = GD.create_county_dist_data(
                GD.county, GD.county_fips, GD.county_dist[cou]
            )
        self.statusbar.showMessage("Figure generated.", TIMEOUT_MS)
        self.fig_distance_gds_plot = GD.create_fig(data)

    # endregion

    # region page_distance_gds_plot
    @Slot()
    def on_pushButton_distance_gds_plot_clicked(self):
        self.fig_distance_gds_plot.show()

    # endregion

    # region page_distance_ind
    @Slot()
    def on_pushButton_distance_ind_clicked(self):
        self.statusbar.showMessage("Validating input...")
        if fieldCheckValid(self.lineEdit_distance_ind):
            self.statusbar.showMessage("Calculating...")
            self.statusbar.showMessage("Figure generated.", TIMEOUT_MS)
            if os.path.isfile(self.lineEdit_distance_ind.text()):
                ind_line = GD.indInput(self.lineEdit_distance_ind.text())
                ind_dist = GD.indInves(GD.county, GD.county_fips, GD.alFreq, ind_line)
                self.fig_distance_ind_plot = GD.create_fig(ind_dist)
                self.stackedWidget_distance.setCurrentIndex(2)
            else:
                self.statusbar.showMessage(
                    "Invalid file name: " + self.lineEdit_distance_ind.text(),
                    TIMEOUT_MS,
                )
        else:
            self.statusbar.showMessage(
                "Invalid file name: " + self.lineEdit_distance_ind.text(), TIMEOUT_MS
            )

    # endregion

    # region page_distance_ind_plot
    @Slot()
    def on_pushButton_distance_ind_plot_clicked(self):
        self.fig_distance_ind_plot.show()

    # endregion

    # endregion

    def closeEvent(self, event):
        global log_window
        if len(opened_windows) == 1:
            ans = QMessageBox.question(
                self,
                "Quit",
                "Would you like to exit this program?",
                QMessageBox.Yes,
                QMessageBox.No,
            )
            if ans == QMessageBox.Yes:
                opened_windows.remove(self)
                if not log_window is None:
                    assert isinstance(log_window, LogWindow)
                    log_window.close()
                event.accept()
            else:
                event.ignore()
        else:
            opened_windows.remove(self)
            event.accept()

    def formCheckValid(self, form: QFormLayout):
        for i in range(form.rowCount()):
            field = form.itemAt(i, QFormLayout.FieldRole).widget()
            if not fieldCheckValid(field):
                self.statusbar.showMessage(
                    "Invalid field: %s"
                    % form.itemAt(i, QFormLayout.LabelRole).widget().text()
                )
                return False
        return True

    def clearPageUpdateUpdate(self):
        if not self.widget_update_update is None:
            self.widget_update_update.setParent(None)
            self.widget_update_update.deleteLater()
            self.widget_update_update = None
        if not self.pushButton_update_update is None:
            self.pushButton_update_update.setParent(None)
            self.pushButton_update_update.deleteLater()
            self.pushButton_update_update = None

    def clearSearchResultSnpFigure(self):
        if not self.figure_search_result is None:
            self.figure_search_result.setParent(None)
            self.figure_search_result.deleteLater()
            self.figure_search_result.close()
            self.figure_search_result = None


class LogWindow(QWidget):
    def __init__(self, parent=None):
        global log_window
        super().__init__(parent)
        self.setupUI()
        log_window = self

    def setupUI(self):
        self.resize(400, 300)
        self.setWindowTitle("Log")

        self.VLayout = QVBoxLayout(self)
        self.setLayout(self.VLayout)

        self.textBrowser = QTextBrowser(self)
        self.textBrowser.setStyleSheet('font: 10pt "Consolas";')
        self.VLayout.addWidget(self.textBrowser)

        self.HLayoutWidget = QWidget(self)
        self.HLayout = QHBoxLayout(self.HLayoutWidget)
        self.VLayout.addWidget(self.HLayoutWidget)

        self.pushButton_refresh = QPushButton(self)
        self.pushButton_refresh.setObjectName("pushButton_refresh")
        self.HLayout.addWidget(self.pushButton_refresh)
        self.pushButton_refresh.setText("Refresh (Ctrl+R)")

        self.pushButton_export = QPushButton(self.HLayoutWidget)
        self.pushButton_export.setObjectName("pushButton_export")
        self.HLayout.addWidget(self.pushButton_export)
        self.pushButton_export.setText("Export")

        self.pushButton_close = QPushButton(self.HLayoutWidget)
        self.pushButton_close.setObjectName("pushButton_close")
        self.HLayout.addWidget(self.pushButton_close)
        self.pushButton_close.setText("Close (Esc)")

        QMetaObject.connectSlotsByName(self)
        self.pushButton_refresh.setShortcut("Ctrl+R")
        self.pushButton_close.setShortcut("Esc")
        self.on_pushButton_refresh_clicked()

    @Slot()
    def on_pushButton_refresh_clicked(self):
        self.textBrowser.setText(logToStr())

    @Slot()
    def on_pushButton_export_clicked(self):
        exportLog()

    @Slot()
    def on_pushButton_close_clicked(self):
        self.close()

    def closeEvent(self, event):
        global log_window
        assert log_window == self
        log_window = None
        event.accept()


class SQLwindow(QDialog):
    def __init__(self):
        super().__init__()
        self.setupUI()

    def setupUI(self):
        self.resize(400, 300)
        self.setWindowTitle("SQL query")

        self.VLayout = QVBoxLayout(self)
        self.setLayout(self.VLayout)

        self.plainTextEdit_SQL = QPlainTextEdit(self)
        self.plainTextEdit_SQL.setStyleSheet('font: 10pt "Consolas";')
        self.VLayout.addWidget(self.plainTextEdit_SQL)

        self.pushButton_SQL = QPushButton(self)
        self.pushButton_SQL.setObjectName("pushButton_SQL")
        self.VLayout.addWidget(self.pushButton_SQL)
        self.pushButton_SQL.setText("Submit")

        QMetaObject.connectSlotsByName(self)

    @Slot()
    def on_pushButton_SQL_clicked(self):
        result = getResponse(self.plainTextEdit_SQL.toPlainText())
        self.plainTextEdit_SQL.appendPlainText(str(result))


print("Launching BioJoin...")
app = QApplication()
MainWindow().show()
app.exec()
print("Program terminated.")
