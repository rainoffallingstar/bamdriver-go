package bamnative

import (
	"testing"
)

func TestSAMFlagConstants(t *testing.T) {
	flags := map[string]uint16{
		"FlagPaired":      FlagPaired,
		"FlagProperPair":   FlagProperPair,
		"FlagUnmapped":     FlagUnmapped,
		"FlagMateUnmapped": FlagMateUnmapped,
	}

	seen := make(map[uint16]bool)
	for name, f := range flags {
		if f == 0 {
			t.Errorf("%s is zero", name)
		}
		if seen[f] {
			t.Errorf("duplicate flag value %#x for %s", f, name)
		}
		seen[f] = true
		if f&(f-1) != 0 {
			t.Errorf("%s = %#x is not a single bit", name, f)
		}
	}
}

func TestRecord_IsPaired(t *testing.T) {
	var rec Record

	rec.Flags = FlagUnmapped
	if rec.IsPaired() {
		t.Error("unmapped record should not be paired")
	}

	rec.Flags = FlagPaired
	if !rec.IsPaired() {
		t.Error("paired record should be paired")
	}

	rec.Flags = FlagPaired | FlagProperPair
	if !rec.IsPaired() {
		t.Error("properly paired record should be paired")
	}
}

func TestRecord_IsUnmapped(t *testing.T) {
	var rec Record

	rec.Flags = 0
	if rec.IsUnmapped() {
		t.Error("mapped record should not be unmapped")
	}

	rec.Flags = FlagUnmapped
	if !rec.IsUnmapped() {
		t.Error("unmapped record should be unmapped")
	}
}

func TestRecord_IsMateUnmapped(t *testing.T) {
	var rec Record

	rec.Flags = 0
	if rec.IsMateUnmapped() {
		t.Error("mate-mapped record should not have mate unmapped")
	}

	rec.Flags = FlagMateUnmapped
	if !rec.IsMateUnmapped() {
		t.Error("mate-unmapped record should have mate unmapped")
	}
}

func TestRecord_GetAuxField(t *testing.T) {
	var rec Record
	aux := rec.GetAuxField("NM")
	if aux != nil {
		t.Log("Unexpected aux field on empty record")
	}
}

func TestAuxField_Type(t *testing.T) {
	var af AuxField
	if af.Tag != "" {
		t.Log("empty AuxField should have empty tag")
	}
}
