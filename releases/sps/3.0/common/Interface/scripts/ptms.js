/**
 * Form field handler for PTM type selection radio buttons.
 */
var PTMTypeFormFieldHandler = {
	getFieldValue: radioFormFieldHandler.getFieldValue,
	setFieldValue: radioFormFieldHandler.setFieldValue,
	
	clearFieldValue: function(field) {
		if (field == null)
			return;
		else if (field.value == "opt")
			field.checked = true;
		else field.checked = false;
	}
}

/**
 * Form field handler for custom PTM parameters.
 */
var customPTMFormFieldHandler = {
	getFieldValue: defaultFormFieldHandler.getFieldValue,
	
	setFieldValue: function(field, value) {
		if (field == null || value == null)
			return;
		else {
			var customPTMs = getFields(field.form, field.name);
			for (var i in customPTMs)
				if (customPTMs[i].value == value)
					return;
			var values = value.split(",");
			if (values != null && values.length == 3) {
				setFieldValue(field.form, "ptm_mass", values[0]);
				setFieldValue(field.form, "ptm_residue", values[1]);
				setFieldValue(field.form, "ptm_type", values[2]);
				addPTM(field.form);
			}
		}
	},
	
	clearFieldValue: function(field) {
		if (field == null)
			return;
		else field.parentNode.cells[0].firstChild.onclick();
	}
}

// assign PTM-related form field handlers
formFieldHandlers["ptm_type"] = PTMTypeFormFieldHandler;
formFieldHandlers["ptm.custom_PTM"] = customPTMFormFieldHandler;

var PTM_DISPLAY = {
	fix: "FIXED",
	opt: "OPTIONAL",
	cterminal: "C-TERMINAL",
	nterminal: "N-TERMINAL"
}
var floatRegExp = /^[\+\-]?(\d+(\.\d*)?|\.?\d+)$/;
var aminoAcidRegExp = /^(\*|[ACDEFGHIKLMNPQRSTVWY]+)$/;

function addPTM(form) {
	var mass = getFieldValue(form, "ptm_mass");
	var resd = getFieldValue(form, "ptm_residue");
	var type = getFieldValue(form, "ptm_type");
	if (!floatRegExp.test(mass)) {
		alert("Mass must be a real number.");
		form.ptm_mass.focus();
		return;
	}
	if (!aminoAcidRegExp.test(resd)) {
		alert("Residues must be an asterisk or " +
			"a string of amino acid abbreviations.");
		form.ptm_residue.focus();
		return;
	}
	var tbl = document.getElementById("PTM_table");
	var row = tbl.insertRow(tbl.rows.length - 1);
	var ctrlCell = row.insertCell(0);
	var massCell = row.insertCell(1);
	var resdCell = row.insertCell(2);
	var typeCell = row.insertCell(3);
	var ctrlButton = document.createElement("img");
	ctrlButton.src = "images/minus.png";
	ctrlButton.className = "selectable"
	ctrlButton.onclick = function() {
		tbl.deleteRow(row.rowIndex);
		return false;
	}
	ctrlCell.appendChild(ctrlButton);
	massCell.innerHTML = mass;
	resdCell.innerHTML = resd;
	typeCell.innerHTML = PTM_DISPLAY[type];
	var param = document.createElement("input");
	param.type = "hidden";
	param.name = "Default.ptm.custom_PTM";
	param.value = mass + "," + resd + "," + type;
	row.appendChild(param);
	clearFieldValue(form, "ptm_mass");
	clearFieldValue(form, "ptm_residue");
	clearFieldValue(form, "ptm_type");
	form.ptm_mass.focus();
}
